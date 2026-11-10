#!/usr/bin/env python3
"""Machine-check scRNA/CM02 production and validation outputs."""

from __future__ import annotations

import os
import random
import sys
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

MASTER_SEED = int(os.environ.get("LIUF_MASTER_SEED", "20260807"))


def initialize_rng() -> None:
    if os.environ.get("PYTHONHASHSEED") != str(MASTER_SEED):
        raise RuntimeError("PYTHONHASHSEED is not the audited master seed.")
    random.seed(MASTER_SEED)
    np.random.seed(MASTER_SEED)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_full_count_audit(
    audit_root: Path,
    result_root: Path,
) -> tuple[pd.DataFrame, dict[str, object]]:
    raw_qc_path = audit_root / "scrna_raw_qc_counts.tsv"
    atlas_metadata_path = (
        result_root
        / "atlas_panel_data"
        / "figure2_atlas_umap_metadata.tsv.gz"
    )
    raw_qc = pd.read_csv(raw_qc_path, sep="\t")
    atlas = pd.read_csv(
        atlas_metadata_path,
        sep="\t",
        usecols=[
            "cell",
            "sample",
            "group",
            "source",
            "celltype",
            "sub_celltype",
        ],
    )
    atlas[["sample", "group", "source", "celltype", "sub_celltype"]] = (
        atlas[
            ["sample", "group", "source", "celltype", "sub_celltype"]
        ].fillna("unassigned")
    )

    prohibited = {"CC1", "CC2", "CC3"}
    prohibited_hits = sorted(prohibited.intersection(set(atlas["sample"])))
    if prohibited_hits:
        raise RuntimeError(
            "Final raw-derived atlas contains prohibited samples: "
            + ",".join(prohibited_hits)
        )
    if atlas["cell"].duplicated().any():
        raise RuntimeError("Final raw-derived atlas metadata has duplicate cells.")

    raw_samples = set(raw_qc["sample"])
    final_samples = set(atlas["sample"])
    if raw_samples != final_samples:
        raise RuntimeError(
            "Raw-QC and final-atlas sample sets differ: "
            f"raw_only={sorted(raw_samples - final_samples)};"
            f"final_only={sorted(final_samples - raw_samples)}"
        )

    final_by_sample = (
        atlas.groupby(["sample", "group", "source"], as_index=False)
        .size()
        .rename(columns={"size": "final_atlas_cells"})
    )
    sample_counts = raw_qc.merge(
        final_by_sample,
        on=["sample", "group", "source"],
        how="left",
        validate="one_to_one",
    )
    sample_counts["final_atlas_cells"] = (
        sample_counts["final_atlas_cells"].fillna(0).astype(int)
    )
    if (
        sample_counts["final_atlas_cells"]
        > sample_counts["post_qc_cells"]
    ).any():
        bad = sample_counts.loc[
            sample_counts["final_atlas_cells"]
            > sample_counts["post_qc_cells"],
            "sample",
        ].tolist()
        raise RuntimeError(
            "Final atlas exceeds post-QC cells for samples: " + ",".join(bad)
        )

    records: list[dict[str, object]] = []

    def append_count(
        stage: str,
        aggregation: str,
        n_cells: int,
        denominator_cells: int,
        provenance: str,
        sample: str = "ALL",
        group: str = "ALL",
        source: str = "ALL",
        label: str = "ALL",
    ) -> None:
        records.append(
            {
                "stage": stage,
                "aggregation": aggregation,
                "sample": sample,
                "group": group,
                "source": source,
                "label": label,
                "n_cells": int(n_cells),
                "denominator_cells": int(denominator_cells),
                "retention_fraction": (
                    float(n_cells) / float(denominator_cells)
                    if denominator_cells > 0
                    else np.nan
                ),
                "input_provenance": provenance,
                "status": "PASS",
            }
        )

    raw_total = int(sample_counts["pre_qc_cells"].sum())
    post_qc_total = int(sample_counts["post_qc_cells"].sum())
    final_total = int(sample_counts["final_atlas_cells"].sum())
    append_count(
        "raw_input",
        "total",
        raw_total,
        raw_total,
        "scrna_raw_qc_counts.tsv:pre_qc_cells",
    )
    append_count(
        "post_qc",
        "total",
        post_qc_total,
        raw_total,
        "scrna_raw_qc_counts.tsv:post_qc_cells",
    )
    append_count(
        "final_atlas",
        "total",
        final_total,
        post_qc_total,
        "figure2_atlas_umap_metadata.tsv.gz",
    )
    for row in sample_counts.itertuples(index=False):
        append_count(
            "raw_input",
            "sample",
            row.pre_qc_cells,
            row.pre_qc_cells,
            "scrna_raw_qc_counts.tsv:pre_qc_cells",
            sample=row.sample,
            group=row.group,
            source=row.source,
        )
        append_count(
            "post_qc",
            "sample",
            row.post_qc_cells,
            row.pre_qc_cells,
            "scrna_raw_qc_counts.tsv:post_qc_cells",
            sample=row.sample,
            group=row.group,
            source=row.source,
        )
        append_count(
            "final_atlas",
            "sample",
            row.final_atlas_cells,
            row.post_qc_cells,
            "figure2_atlas_umap_metadata.tsv.gz",
            sample=row.sample,
            group=row.group,
            source=row.source,
        )

    for aggregation, columns in (
        ("group", ["group"]),
        ("celltype", ["group", "celltype"]),
        ("sub_celltype", ["group", "celltype", "sub_celltype"]),
    ):
        grouped = atlas.groupby(columns, as_index=False).size()
        for row in grouped.itertuples(index=False):
            values = dict(zip(columns, row[:-1], strict=True))
            label_parts = [
                str(values[key])
                for key in ("celltype", "sub_celltype")
                if key in values
            ]
            append_count(
                "final_atlas",
                aggregation,
                row[-1],
                final_total,
                "figure2_atlas_umap_metadata.tsv.gz",
                group=str(values.get("group", "ALL")),
                label=" / ".join(label_parts) if label_parts else "ALL",
            )

    count_audit = pd.DataFrame(records)
    count_audit.to_csv(
        audit_root / "scrna_full_count_rebuild.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    metrics: dict[str, object] = {
        "raw_total": raw_total,
        "post_qc_total": post_qc_total,
        "final_total": final_total,
        "n_samples": len(final_samples),
        "prohibited_hits": len(prohibited_hits),
        "audit_rows": len(count_audit),
    }
    return count_audit, metrics


def main() -> None:
    initialize_rng()
    mode = sys.argv[1] if len(sys.argv) > 1 else "audit"
    if mode not in {"audit", "curated", "full"}:
        raise ValueError("Mode must be one of: audit, curated, full")
    root = Path(__file__).resolve().parents[3]
    audit_root = root / "audit"
    result_root = root / "results" / "20_scrna_cm02"
    pipeline_root = root / "code" / "pipelines" / "20_scrna_cm02"
    rows: list[dict[str, object]] = []

    strict_seed_path = audit_root / "scrna_seed_audit.tsv"
    seed_plan_path = audit_root / "scrna_seed_audit_publication_plan.tsv"
    expected_columns = {
        "pipeline",
        "step",
        "script_path",
        "language",
        "algorithm",
        "seed",
        "seed_source",
        "rng_kind",
        "status",
        "notes",
    }
    plan = pd.read_csv(seed_plan_path, sep="\t")
    plan_ok = {
        "pipeline",
        "step",
        "script_path",
        "language",
        "algorithm",
        "seed",
        "seed_source",
        "rng_kind",
        "planned_status",
        "publication_gate",
        "notes",
    }.issubset(plan.columns)
    rows.append(
        {
            "check": "seed publication plan available",
            "observed": len(plan),
            "expected": "planned rows with publication gates",
            "status": "PASS" if plan_ok else "FAIL",
        }
    )
    repeat_path = audit_root / "scrna_seed_repeatability.tsv"
    strict_withheld_ok = (
        mode in {"curated", "full"} or not strict_seed_path.exists()
    )
    rows.append(
        {
            "check": "raw-atlas seed audit policy",
            "observed": strict_seed_path.exists(),
            "expected": (
                "not required by frozen-curated production"
                if mode == "curated"
                else (
                    "published by completed 09"
                    if mode == "full"
                    else False
                )
            ),
            "status": "PASS" if strict_withheld_ok else "FAIL",
        }
    )

    cc = pd.read_csv(audit_root / "scrna_cc_exclusion.tsv", sep="\t")
    cc_hits = int(cc["prohibited_hits"].sum())
    rows.append(
        {
            "check": "CC1-CC3 excluded",
            "observed": cc_hits,
            "expected": 0,
            "status": "PASS" if cc_hits == 0 else "FAIL",
        }
    )

    curation_audit_path = audit_root / "scrna_frozen_curation_objects.tsv"
    curation = pd.read_csv(curation_audit_path, sep="\t")
    expected_curation_objects = {
        "atlas",
        "hepatocyte",
        "t_cell",
        "macrophage",
        "fibroblast",
        "b_cell",
    }
    curation_ok = (
        set(curation["object"]) == expected_curation_objects
        and curation["status"].eq("PASS").all()
        and curation["prohibited_cc_hits"].eq(0).all()
        and curation["sha256"].eq(curation["expected_sha256"]).all()
    )
    rows.append(
        {
            "check": "frozen expert-curation boundary",
            "observed": ",".join(sorted(curation["object"].astype(str))),
            "expected": "six checksum-locked PASS objects; CC1-CC3=0",
            "status": "PASS" if curation_ok else "FAIL",
        }
    )

    forbidden_tokens: list[str] = []
    for filename in (
        "03_rebuild_atlas_from_local_raw.R",
        "04_rebuild_atlas_panel_data.R",
        "05_rebuild_covarnet_cm02.R",
    ):
        text = (pipeline_root / filename).read_text(encoding="utf-8")
        for token in (
            "validation_root",
            "frozen_targets",
            "combined_annotated_v5.rds",
            "cm_network.rds",
        ):
            if token in text:
                forbidden_tokens.append(f"{filename}:{token}")
    rows.append(
        {
            "check": "production scripts do not read validation targets",
            "observed": ";".join(forbidden_tokens),
            "expected": "",
            "status": "PASS" if not forbidden_tokens else "FAIL",
        }
    )

    if mode in {"curated", "full"}:
        curated_outputs = [
            result_root / "covarnet" / "freq_inputs_frozen_curated.rds",
            result_root / "covarnet" / "NMF_K5_frozen_curated.rds",
            result_root / "covarnet" / "cm_network_frozen_curated.rds",
            audit_root / "scrna_frozen_curation_objects.tsv",
            audit_root / "scrna_frozen_curated_vs_validation_frequency.tsv",
            audit_root / "scrna_frozen_curated_vs_validation_nmf.tsv",
            audit_root / "scrna_frozen_curated_vs_validation_network.tsv",
        ]
        missing_curated = [
            str(path) for path in curated_outputs if not path.exists()
        ]
        rows.append(
            {
                "check": "frozen-curated CoVarNet outputs exist",
                "observed": ";".join(missing_curated),
                "expected": "",
                "status": "PASS" if not missing_curated else "FAIL",
            }
        )
        if not missing_curated:
            frequency = pd.read_csv(
                audit_root
                / "scrna_frozen_curated_vs_validation_frequency.tsv",
                sep="\t",
            )
            network = pd.read_csv(
                audit_root
                / "scrna_frozen_curated_vs_validation_network.tsv",
                sep="\t",
            )
            nmf = pd.read_csv(
                audit_root / "scrna_frozen_curated_vs_validation_nmf.tsv",
                sep="\t",
            )
            comparison_ok = (
                frequency["status"].eq("COMPARISON_COMPLETE").all()
                and network["status"].eq("COMPARISON_COMPLETE").all()
                and not nmf.empty
                and nmf["cosine"].notna().all()
            )
            rows.append(
                {
                    "check": "frozen-curated CoVarNet comparison complete",
                    "observed": (
                        f"frequency_rows={len(frequency)};"
                        f"nmf_pairs={len(nmf)};network_rows={len(network)}"
                    ),
                    "expected": "all comparison receipts complete",
                    "status": "PASS" if comparison_ok else "FAIL",
                }
            )

    if mode == "full":
        count_audit, count_metrics = build_full_count_audit(
            audit_root,
            result_root,
        )
        rows.append(
            {
                "check": "full raw-derived count/QC audit",
                "observed": (
                    f"raw={count_metrics['raw_total']};"
                    f"post_qc={count_metrics['post_qc_total']};"
                    f"final={count_metrics['final_total']};"
                    f"samples={count_metrics['n_samples']};"
                    f"rows={count_metrics['audit_rows']}"
                ),
                "expected": "39 samples; CC1-CC3=0; final<=post_qc",
                "status": (
                    "PASS"
                    if count_metrics["n_samples"] == 39
                    and count_metrics["prohibited_hits"] == 0
                    and not count_audit.empty
                    else "FAIL"
                ),
            }
        )
        if not strict_seed_path.exists():
            raise RuntimeError(
                "Strict seed audit is absent; repeatability publication "
                "has not completed."
            )
        seed = pd.read_csv(strict_seed_path, sep="\t")
        seed_ok = expected_columns.issubset(seed.columns) and not seed[
            list(expected_columns)
        ].isna().any().any() and set(seed["status"]).issubset(
            {"fixed", "deterministic", "blocked"}
        )
        rows.append(
            {
                "check": "strict seed audit complete",
                "observed": len(seed),
                "expected": "all fields populated",
                "status": "PASS" if seed_ok else "FAIL",
            }
        )

        dependency_path = audit_root / "scrna_seed_audit_dependencies.tsv"
        dependency = pd.read_csv(dependency_path, sep="\t")
        dependency_records_ok = True
        dependency_errors: list[str] = []
        seed_mtime = strict_seed_path.stat().st_mtime
        for record in dependency.itertuples(index=False):
            source_path = root / record.path
            if not source_path.exists():
                dependency_records_ok = False
                dependency_errors.append(f"missing:{record.path}")
                continue
            observed_hash = sha256_file(source_path)
            if observed_hash != record.sha256:
                dependency_records_ok = False
                dependency_errors.append(f"hash:{record.path}")
            if seed_mtime <= source_path.stat().st_mtime:
                dependency_records_ok = False
                dependency_errors.append(f"mtime:{record.path}")
        rows.append(
            {
                "check": "strict audit regenerated from current dependencies",
                "observed": ";".join(dependency_errors),
                "expected": "",
                "status": "PASS" if dependency_records_ok else "FAIL",
            }
        )

        repeat = pd.read_csv(repeat_path, sep="\t")
        repeat_complete = {
            "FindNeighbors/Annoy SNN",
            "FindClusters/Louvain",
            "RunUMAP/uwot",
            "NMF/nsNMF K=5",
        }.issubset(set(repeat["algorithm"]))
        rows.append(
            {
                "check": "Annoy Louvain UMAP NMF duplicate tests completed",
                "observed": ",".join(sorted(set(repeat["algorithm"]))),
                "expected": "all four algorithm classes",
                "status": "PASS" if repeat_complete else "FAIL",
            }
        )
        expected_outputs = [
            result_root / "atlas_major_raw_derived.rds",
            result_root / "atlas_final_raw_derived.rds",
            result_root / "hep_raw_derived.rds",
            result_root / "t_raw_derived.rds",
            result_root / "mp_raw_derived.rds",
            result_root / "fib_raw_derived.rds",
            result_root / "b_raw_derived.rds",
            result_root / "covarnet" / "NMF_K5_frozen_curated.rds",
            result_root / "covarnet" / "cm_network_frozen_curated.rds",
            audit_root / "scrna_raw_vs_frozen_objects.tsv",
            audit_root / "scrna_frozen_curated_vs_validation_nmf.tsv",
            audit_root / "scrna_frozen_curation_objects.tsv",
            audit_root / "scrna_full_count_rebuild.tsv",
        ]
        missing = [str(path) for path in expected_outputs if not path.exists()]
        rows.append(
            {
                "check": "raw-derived audit and frozen-curated outputs exist",
                "observed": ";".join(missing),
                "expected": "",
                "status": "PASS" if not missing else "FAIL",
            }
        )
        composition = pd.read_csv(
            result_root
            / "atlas_panel_data"
            / "figure2_composition_values.tsv",
            sep="\t",
        )
        sums = composition.groupby("group")[
            "proportion_within_group"
        ].sum()
        maximum_error = float(np.max(np.abs(sums.to_numpy() - 1.0)))
        rows.append(
            {
                "check": "raw-derived Figure2 composition audit sums",
                "observed": maximum_error,
                "expected": "<=1e-12",
                "status": "PASS" if maximum_error <= 1e-12 else "FAIL",
            }
        )

    validation = pd.DataFrame(rows)
    validation.to_csv(
        audit_root / "scrna_output_validation.tsv",
        sep="\t",
        index=False,
    )
    print(validation.to_string(index=False))
    if validation["status"].ne("PASS").any():
        raise RuntimeError("scRNA/CM02 output validation failed.")


if __name__ == "__main__":
    main()
