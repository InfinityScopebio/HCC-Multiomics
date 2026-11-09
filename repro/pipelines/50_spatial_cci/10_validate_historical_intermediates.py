#!/usr/bin/env python3
"""Validate checksum-locked historical intermediates used downstream."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from config import load_config


REQUIRED_COLUMNS = (
    "role",
    "section",
    "source_path",
    "materialized_path",
    "bytes",
    "sha256",
    "article_use",
)


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_schema(table: pd.DataFrame, sections: tuple[str, ...]) -> None:
    if tuple(table.columns) != REQUIRED_COLUMNS:
        raise RuntimeError(
            "Historical-intermediate manifest schema mismatch: "
            f"{tuple(table.columns)}"
        )
    if table[["role", "section"]].duplicated().any():
        raise RuntimeError("Historical-intermediate manifest has duplicate keys.")
    observed_sections = set(
        table.loc[table["role"].eq("rctd_weight"), "section"].astype(str)
    )
    if observed_sections != set(sections):
        raise RuntimeError(
            "Historical RCTD section set mismatch; missing="
            f"{sorted(set(sections) - observed_sections)}, unexpected="
            f"{sorted(observed_sections - set(sections))}"
        )
    singleton_roles = {"rctd_reference", "covarnet_model"}
    for role in singleton_roles:
        rows = table[table["role"].eq(role)]
        if len(rows) != 1 or rows.iloc[0]["section"] != "ALL":
            raise RuntimeError(f"Manifest requires one ALL row for {role}.")
    if set(table["role"]) != {"rctd_reference", "rctd_weight", "covarnet_model"}:
        raise RuntimeError("Historical-intermediate manifest has an unknown role.")


def main() -> None:
    config = load_config()
    manifest_path = (
        config.project_root
        / "repro/manifests/spatial_historical_intermediates.tsv"
    )
    table = pd.read_csv(manifest_path, sep="\t", dtype={"section": str})
    validate_schema(table, config.all_sections)

    rows: list[dict[str, object]] = []
    for record in table.to_dict(orient="records"):
        path = config.assert_production_input(
            config.project_root / str(record["materialized_path"])
        )
        exists = path.is_file()
        observed_bytes = path.stat().st_size if exists else -1
        observed_sha256 = sha256_file(path) if exists else ""
        expected_bytes = int(record["bytes"])
        expected_sha256 = str(record["sha256"])
        status = (
            "PASS"
            if exists
            and observed_bytes == expected_bytes
            and observed_sha256 == expected_sha256
            else "FAIL"
        )
        rows.append(
            {
                **record,
                "resolved_path": str(path),
                "exists": exists,
                "observed_bytes": observed_bytes,
                "observed_sha256": observed_sha256,
                "status": status,
            }
        )

    audit = pd.DataFrame(rows)
    output = (
        config.audit_root
        / "spatial_historical_intermediate_validation.tsv"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output, sep="\t", index=False)
    failures = audit[~audit["status"].eq("PASS")]
    if not failures.empty:
        raise RuntimeError(
            "Historical intermediate validation failed for: "
            + ", ".join(failures["materialized_path"].astype(str))
        )
    print(
        "Historical article-equivalence intermediates ready: "
        f"{len(audit)} files, {len(config.all_sections)} RCTD sections."
    )


if __name__ == "__main__":
    main()
