#!/usr/bin/env python3
"""Audit local raw scRNA inputs without reading validation targets."""

from __future__ import annotations

import hashlib
import os
import random
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

MASTER_SEED = int(os.environ.get("LIUF_MASTER_SEED", "20260807"))
ELIGIBLE_INTERNAL = ("HCC1", "HCC2", "HCC3", "NC1", "NC2", "NC3")
PROHIBITED_INTERNAL = ("CC1", "CC2", "CC3")


@dataclass(frozen=True)
class Paths:
    project_root: Path
    result_root: Path
    audit_root: Path
    internal_root: Path
    gse115469: Path
    gse136103: Path
    gse149614_counts: Path
    gse149614_metadata: Path


def resolve_paths() -> Paths:
    pipeline_dir = Path(__file__).resolve().parent
    project_root = pipeline_dir.parents[2]
    return Paths(
        project_root=project_root,
        result_root=project_root / "results" / "20_scrna_cm02",
        audit_root=project_root / "audit",
        internal_root=(
            project_root
            / "repo"
            / "data"
            / "private"
            / "raw"
            / "data"
            / "scRNAData"
        ),
        gse115469=(
            project_root
            / "repo"
            / "data"
            / "public"
            / "raw"
            / "GSE115469"
            / "GSE115469_Data.csv.gz"
        ),
        gse136103=(
            project_root
            / "repo"
            / "data"
            / "public"
            / "raw"
            / "GSE136103"
            / "GSE136103_extracted"
        ),
        gse149614_counts=(
            project_root
            / "repo"
            / "data"
            / "public"
            / "raw"
            / "GSE149614"
            / "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
        ),
        gse149614_metadata=(
            project_root
            / "repo"
            / "data"
            / "public"
            / "raw"
            / "GSE149614"
            / "GSE149614_HCC.metadata.updated.txt.gz"
        ),
    )


def initialize_rng() -> None:
    observed_hash_seed = os.environ.get("PYTHONHASHSEED")
    if observed_hash_seed != str(MASTER_SEED):
        raise RuntimeError(
            "PYTHONHASHSEED must be exported as "
            f"{MASTER_SEED}; observed={observed_hash_seed!r}"
        )
    random.seed(MASTER_SEED)
    np.random.seed(MASTER_SEED)


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def internal_files(paths: Paths) -> list[Path]:
    observed = sorted(
        path.name for path in paths.internal_root.iterdir() if path.is_dir()
    )
    if observed != sorted(ELIGIBLE_INTERNAL):
        raise RuntimeError(
            "Internal raw root must contain exactly HCC1-HCC3 and NC1-NC3; "
            f"observed={observed}"
        )
    prohibited_hits = sorted(set(observed) & set(PROHIBITED_INTERNAL))
    if prohibited_hits:
        raise RuntimeError(f"Prohibited CC inputs found: {prohibited_hits}")

    expected = ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
    files: list[Path] = []
    for sample in ELIGIBLE_INTERNAL:
        matrix_dir = (
            paths.internal_root / sample / "filtered_feature_bc_matrix"
        )
        observed_files = sorted(
            path.name for path in matrix_dir.iterdir() if path.is_file()
        )
        if observed_files != sorted(expected):
            raise RuntimeError(
                f"{sample} 10x files differ: observed={observed_files}"
            )
        files.extend(matrix_dir / filename for filename in expected)
    return files


def public_files(paths: Paths) -> list[tuple[Path, str, str]]:
    records: list[tuple[Path, str, str]] = [
        (paths.gse115469, "included", ""),
        (paths.gse149614_counts, "included", ""),
        (paths.gse149614_metadata, "included", ""),
    ]
    suffixes = (
        "_barcodes.tsv.gz",
        "_genes.tsv.gz",
        "_matrix.mtx.gz",
    )
    for path in sorted(paths.gse136103.iterdir()):
        if not path.is_file() or not path.name.endswith(suffixes):
            continue
        excluded = "_blood" in path.name or "_mouse" in path.name
        records.append(
            (
                path,
                "excluded_before_downstream" if excluded else "included",
                "GSE136103 blood or mouse library" if excluded else "",
            )
        )
    return records


def main() -> None:
    initialize_rng()
    paths = resolve_paths()
    paths.result_root.mkdir(parents=True, exist_ok=True)
    paths.audit_root.mkdir(parents=True, exist_ok=True)

    internal = internal_files(paths)
    public = public_files(paths)
    records = [
        (path, "included", "")
        for path in internal
    ] + public
    for path, _, _ in records:
        if not path.exists():
            raise FileNotFoundError(path)

    inventory = pd.DataFrame(
        [
            {
                "path": str(path.relative_to(paths.project_root)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "usage_status": usage_status,
                "reason": reason,
            }
            for path, usage_status, reason in records
        ]
    )
    inventory.to_csv(
        paths.result_root / "local_scrna_input_inventory.tsv",
        sep="\t",
        index=False,
    )

    cc_audit = pd.DataFrame(
        [
            {
                "stage": "raw_input_directories",
                "eligible_samples": ";".join(ELIGIBLE_INTERNAL),
                "prohibited_samples": ";".join(PROHIBITED_INTERNAL),
                "prohibited_hits": 0,
                "status": "PASS",
            }
        ]
    )
    cc_audit.to_csv(
        paths.audit_root / "scrna_cc_exclusion.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        [
            {"component": "Python", "version": sys.version.split()[0]},
            {"component": "pandas", "version": pd.__version__},
            {"component": "numpy", "version": np.__version__},
            {
                "component": "PYTHONHASHSEED",
                "version": os.environ["PYTHONHASHSEED"],
            },
        ]
    ).to_csv(
        paths.audit_root / "scrna_runtime_python.tsv",
        sep="\t",
        index=False,
    )
    print(
        f"PASS: {len(inventory)} local raw files hashed; "
        "CC1-CC3 excluded; no validation target read."
    )


if __name__ == "__main__":
    main()
