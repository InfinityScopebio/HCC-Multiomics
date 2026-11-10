#!/usr/bin/env python3
"""Summarize frozen-curated scRNA production and raw-derived audit status."""

from __future__ import annotations

import os
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MASTER_SEED = int(os.environ.get("LIUF_MASTER_SEED", "20260807"))


def initialize_rng() -> None:
    if os.environ.get("PYTHONHASHSEED") != str(MASTER_SEED):
        raise RuntimeError("PYTHONHASHSEED is not the audited master seed.")
    random.seed(MASTER_SEED)
    np.random.seed(MASTER_SEED)


def main() -> None:
    initialize_rng()
    mode = sys.argv[1] if len(sys.argv) > 1 else "audit"
    if mode not in {"audit", "curated", "full"}:
        raise ValueError("Mode must be one of: audit, curated, full")
    project_root = Path(__file__).resolve().parents[3]
    audit_root = project_root / "audit"
    result_root = project_root / "results" / "20_scrna_cm02"

    if mode in {"curated", "full"}:
        curation = pd.read_csv(
            audit_root / "scrna_frozen_curation_objects.tsv",
            sep="\t",
        )
        atlas_curation = curation.loc[curation["object"] == "atlas"].iloc[0]
        hepatocyte_curation = curation.loc[
            curation["object"] == "hepatocyte"
        ].iloc[0]
        atlas_validation = f"curated atlas SHA256={atlas_curation['sha256']}"
        if mode == "full":
            object_comparison = pd.read_csv(
                audit_root / "scrna_raw_vs_frozen_objects.tsv",
                sep="\t",
            )
            atlas_row = object_comparison.loc[
                object_comparison["object"] == "atlas"
            ].iloc[0]
            atlas_validation += (
                "; independent raw diagnostic "
                f"cell_jaccard={atlas_row['cell_jaccard']:.6f}"
            )
        frequency = pd.read_csv(
            audit_root / "scrna_frozen_curated_vs_validation_frequency.tsv",
            sep="\t",
        ).iloc[0]
        network = pd.read_csv(
            audit_root / "scrna_frozen_curated_vs_validation_network.tsv",
            sep="\t",
        ).iloc[0]
        statuses = [
            {
                "scope": "Figure2-5_frozen_scRNA_atlas",
                "production_input": (
                    "checksum-locked expert-curated historical atlas"
                ),
                "production_output": (
                    "data/curation/scRNA/combined_annotated_v5.rds"
                ),
                "validation_metric": atlas_validation,
                "status": "FROZEN_EXPERT_CURATION_BOUNDARY",
            },
            {
                "scope": "Figure2_subclusters_and_summaries",
                "production_input": (
                    "checksum-locked expert-curated T/macrophage/fibroblast/"
                    "B/hepatocyte historical objects"
                ),
                "production_output": (
                    "data/curation/scRNA/*_annotated.rds"
                ),
                "validation_metric": "per-object SHA256/dimension/metadata audit",
                "status": "FROZEN_EXPERT_CURATION_BOUNDARY",
            },
            {
                "scope": "Figure3_CM02_CoVarNet",
                "production_input": (
                    "frozen-curated final cell-state frequencies; K=5; "
                    "vendored CoVarNet v0.3.0"
                ),
                "production_output": (
                    "NMF_K5_frozen_curated.rds + "
                    "cm_network_frozen_curated.rds"
                ),
                "validation_metric": (
                    f"frequency_r={frequency['pearson_correlation']:.6f};"
                    f"CM02_network_jaccard={network['jaccard']:.6f}"
                ),
                "status": "FROZEN_CURATED_DERIVED_AND_VALIDATED",
            },
            {
                "scope": "Figure3_bulk_trajectory",
                "production_input": (
                    "raw bulk matrices projected with a MuSiC reference "
                    "built from the frozen-curated atlas"
                ),
                "production_output": "raw_to_figure3/music/sc_eset.rds",
                "validation_metric": (
                    "sc_reference_provenance.tsv records atlas SHA256"
                ),
                "status": "FROZEN_CURATED_REFERENCE",
            },
            {
                "scope": "Figure4_scRNA_derived_signature",
                "production_input": (
                    "checksum-locked Hep_C UP/DOWN signature derived from "
                    "the frozen-curated hepatocyte object"
                ),
                "production_output": (
                    "data/curation/resistance/HepC_signatures_main.gmt"
                ),
                "validation_metric": (
                    "frozen hepatocyte SHA256="
                    f"{hepatocyte_curation['sha256']}"
                ),
                "status": "FROZEN_CURATED_DERIVATIVE",
            },
            {
                "scope": "Figure5_scRNA_profiles_and_reference",
                "production_input": (
                    "frozen-curated atlas and hepatocyte object; historical "
                    "RCTD weights remain the article-equivalence boundary"
                ),
                "production_output": (
                    "frozen-atlas LR profiles, Hep_C4 expression background, "
                    "and CoVarNet reference audit"
                ),
                "validation_metric": "atlas and hepatocyte SHA256 exact",
                "status": "FROZEN_CURATED_REFERENCE",
            },
        ]
    else:
        statuses = [
            {
                "scope": "scRNA frozen-curation production plus raw audit",
                "production_input": "curated objects and local raw matrices audited",
                "production_output": "not generated in audit mode",
                "validation_metric": "validation targets inventoried",
                "status": "READY_FOR_FULL_RUN",
            }
        ]

    status_frame = pd.DataFrame(statuses)
    status_frame.to_csv(
        audit_root / "scrna_panel_rebuild_status.tsv",
        sep="\t",
        index=False,
    )

    seed_audit = pd.read_csv(
        audit_root / "scrna_seed_audit_detailed.tsv",
        sep="\t",
    )
    historical = seed_audit.loc[seed_audit["source"] == "historical"]
    report = f"""# scRNA-seq frozen-curation production and raw-derived audit

## Run mode

`{mode}`

## Integrity rules

- Figure 2-5 annotations and downstream single-cell-derived references use six
  checksum-locked expert-curated historical Seurat objects.
- Local raw counts are reconstructed only in the explicit `full` diagnostic
  branch and are never substituted for historical cell-level
  contamination-removal decisions.
- CC1-CC3 are rejected before any downstream operation.
- Frozen CoVarNet and trajectory targets remain validation-only; CoVarNet is
  recomputed from the declared frozen-curated atlas.
- Figure 6 consumes no single-cell object or single-cell-derived signature.
- Master seed is {MASTER_SEED}; {len(historical)} historical algorithm seeds
  are retained with reasons in `scrna_seed_audit.tsv`.
- CoVarNet is vendored from upstream v0.3.0 commit
  `d4c52a4b2ea821dc4b75861c3772e21232d1ee51`.

## Curation boundary

Historical manual exclusion of mixed and contaminating cells was not fully
encoded as a deterministic cell-level rule. The accepted atlas and five
subcluster objects are therefore an explicit expert-curation boundary, with
immutable hashes and object-level audits. The raw-derived RPCA/subclustering
branch remains available to quantify reconstruction drift without redefining
the accepted Figure 2 cell membership or labels.
"""
    (audit_root / "scrna_reproduction_report.md").write_text(
        report,
        encoding="utf-8",
    )
    (result_root / "AUDIT_COMPLETE").write_text(
        f"scRNA/CM02 {mode} audit completed.\n",
        encoding="ascii",
    )
    print(status_frame.to_string(index=False))


if __name__ == "__main__":
    main()
