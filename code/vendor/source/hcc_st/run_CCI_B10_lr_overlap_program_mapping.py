#!/usr/bin/env python3
"""CCI B10 - Map COMMOT CellChat pathways to programs by LR-gene overlap.

This version does not use NicheNet ligand ranking or spot-level program scores.
Each active CellChat LR pair targeting Hep_C4 is matched directly to the four
program gene sets by LR-gene overlap. An LR matching multiple programs has its
interaction score split across those programs to avoid duplicate attribution.
"""
import json
import os
import re
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd


OUT = "hcc_st/CCI_nichenet_commot"
PROGRAM_GENES = "hcc_st/CCI_nichenet/target_sets_4pathway_expressed_long.csv"
LR_AWARE_FLOW_META = f"{OUT}/true_pairwise_lr_edges.csv"

PROGRAM_ORDER = ["MAPK", "NFkB", "mTOR", "Hypoxia"]
SELECT_TOP_N = 10
MECHANISM_PATHWAYS = {
    "MAPK": {"FGF", "GDF", "PDGF"},
    "NFkB": {"CXCL", "CCL"},
    "mTOR": {"CXCL", "EGF"},
    "Hypoxia": {"ANGPTL", "MK"},
}
PROGRAM_LABELS = {
    "MAPK": "MAPK-ERK",
    "NFkB": "NF-kB",
    "mTOR": "PI3K-mTOR",
    "Hypoxia": "Hypoxia",
}


def split_lr_genes(value: str) -> set[str]:
    """Split CellChat ligand/receptor complexes into gene symbols."""
    tokens = re.split(r"[_+&|:]", str(value))
    return {token.strip().upper() for token in tokens if token.strip()}


def collapse_gene_lists(values: pd.Series) -> str:
    genes = {
        gene
        for value in values
        for gene in str(value).split(";")
        if gene
    }
    return ";".join(sorted(genes))


def load_c4_lr_flow() -> pd.DataFrame:
    if not os.path.exists(LR_AWARE_FLOW_META):
        raise FileNotFoundError(
            f"Missing {LR_AWARE_FLOW_META}; run run_CCI_B11_true_pairwise_networks.py first"
        )
    flow = pd.read_csv(LR_AWARE_FLOW_META)
    required = {
        "ligand",
        "receptor",
        "pathway",
        "sender",
        "sender_label",
        "receiver",
        "mean_patient_flow_share",
        "projected_flow_sum",
        "support_fraction",
    }
    missing = required.difference(flow.columns)
    if missing:
        raise KeyError(f"LR-aware flow table lacks columns: {sorted(missing)}")
    flow = flow[
        (flow["receiver"] == "Hep_C4") & (flow["sender"] != "Hep_C4")
    ].copy()
    flow = flow.reset_index(drop=True)
    flow["flow_row_id"] = np.arange(len(flow), dtype=int)
    return flow


def load_program_genes() -> tuple[dict[str, set[str]], dict[str, dict]]:
    pg = pd.read_csv(PROGRAM_GENES)
    pg["program"] = pg["program"].astype(str)
    pg["gene"] = pg["gene"].astype(str).str.upper()
    program_genes = {program: set(g["gene"]) for program, g in pg.groupby("program")}
    gene_set_meta = (
        pg.groupby("program")
        .agg(
            gene_set_name=("gene_set_name", lambda x: ";".join(sorted(set(map(str, x))))),
            target_filter=("target_filter", lambda x: ";".join(sorted(set(map(str, x))))),
            expressed_target_n=("gene", "nunique"),
        )
        .to_dict("index")
    )
    missing_programs = set(PROGRAM_ORDER).difference(program_genes)
    if missing_programs:
        raise KeyError(f"Program gene table lacks: {sorted(missing_programs)}")
    return program_genes, gene_set_meta


def assign_program_weights(mapping: pd.DataFrame) -> pd.DataFrame:
    weighted = mapping.copy()
    match_count = weighted.groupby("flow_row_id")["program"].nunique()
    weighted["program_match_count"] = weighted["flow_row_id"].map(match_count)
    weighted["assignment_fraction"] = 1.0 / weighted["program_match_count"]
    weighted["pre_split_cci_weight"] = weighted["mean_patient_flow_share"]
    weighted["cci_weight"] = (
        weighted["mean_patient_flow_share"] * weighted["assignment_fraction"]
    )
    weighted["assigned_projected_flow"] = (
        weighted["projected_flow_sum"] * weighted["assignment_fraction"]
    )
    return weighted


def build_active_lr_mapping(flow: pd.DataFrame) -> pd.DataFrame:
    program_genes, gene_set_meta = load_program_genes()
    rows = []
    for lr_row in flow.itertuples(index=False):
        ligand_genes = split_lr_genes(lr_row.ligand)
        receptor_genes = split_lr_genes(lr_row.receptor)
        lr_genes = ligand_genes | receptor_genes
        for program in PROGRAM_ORDER:
            if lr_row.pathway not in MECHANISM_PATHWAYS[program]:
                continue
            genes = program_genes[program]
            overlap = sorted(lr_genes & genes)
            if not overlap:
                continue
            meta = gene_set_meta.get(program, {})
            overlap_role = []
            if ligand_genes & genes:
                overlap_role.append("ligand")
            if receptor_genes & genes:
                overlap_role.append("receptor")
            record = lr_row._asdict()
            record.update(
                {
                    "program": program,
                    "program_label": PROGRAM_LABELS[program],
                    "lr_pair": f"{lr_row.ligand}->{lr_row.receptor}",
                    "n_lr_genes": len(lr_genes),
                    "n_program_genes": len(genes),
                    "overlap_n": len(overlap),
                    "overlap_frac_lr": len(overlap) / len(lr_genes),
                    "overlap_frac_program": len(overlap) / len(genes),
                    "jaccard": len(overlap) / len(lr_genes | genes),
                    "overlap_genes": ";".join(overlap),
                    "overlap_role": "+".join(overlap_role),
                    "lr_genes": ";".join(sorted(lr_genes)),
                    "program_gene_set": meta.get("gene_set_name", ""),
                    "target_filter": meta.get("target_filter", ""),
                    "expressed_target_n": meta.get("expressed_target_n", len(genes)),
                    "gene_set_source": os.path.basename(PROGRAM_GENES),
                }
            )
            rows.append(record)

    mapping = pd.DataFrame(rows)
    if mapping.empty:
        raise RuntimeError("No active Hep_C4-directed LR pair overlaps the mechanism programs")
    mapping = assign_program_weights(mapping)
    mapping["flow_weight_source"] = "active_lr_expression_gated_mean_patient_score_share"
    mapping["mapping_source"] = "active_cellchat_lr_gene_overlap_expressed_targets_no_DE"
    mapping["program_order"] = mapping["program"].map(
        {program: index for index, program in enumerate(PROGRAM_ORDER)}
    )
    mapping = mapping.sort_values(
        ["program_order", "cci_weight", "support_fraction", "overlap_n"],
        ascending=[True, False, False, False],
    )
    mapping.to_csv(f"{OUT}/lr_overlap_active_lr_program_mapping.csv", index=False)

    conservation = (
        mapping.groupby("flow_row_id", as_index=False)
        .agg(
            ligand=("ligand", "first"),
            receptor=("receptor", "first"),
            pathway=("pathway", "first"),
            sender=("sender", "first"),
            matched_programs=("program", lambda values: ";".join(sorted(values))),
            program_match_count=("program", "nunique"),
            assignment_fraction_sum=("assignment_fraction", "sum"),
            pre_split_cci_weight=("pre_split_cci_weight", "first"),
            assigned_cci_weight_sum=("cci_weight", "sum"),
        )
    )
    conservation["weight_conserved"] = np.isclose(
        conservation["pre_split_cci_weight"],
        conservation["assigned_cci_weight_sum"],
    )
    conservation.to_csv(f"{OUT}/lr_overlap_assignment_conservation.csv", index=False)

    pathway_mapping = (
        mapping.groupby(["program", "program_label", "pathway"], as_index=False)
        .agg(
            n_active_lr_rows=("flow_row_id", "nunique"),
            n_senders=("sender", "nunique"),
            assigned_cci_weight=("cci_weight", "sum"),
            overlap_genes=("overlap_genes", collapse_gene_lists),
        )
        .sort_values(["program", "assigned_cci_weight"], ascending=[True, False])
    )
    pathway_mapping.to_csv(f"{OUT}/lr_overlap_program_pathway_mapping.csv", index=False)
    pathway_mapping.to_csv(
        f"{OUT}/lr_overlap_program_pathway_mapping_nonzero.csv", index=False
    )
    return mapping


def select_c4_targets(mapping: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = mapping.copy()
    candidates = candidates.sort_values(
        [
            "program_order",
            "cci_weight",
            "support_fraction",
            "overlap_n",
            "jaccard",
        ],
        ascending=[True, False, False, False, False],
    )
    candidates["rank_within_program"] = candidates.groupby("program").cumcount() + 1
    candidates["selected"] = candidates["rank_within_program"] <= SELECT_TOP_N
    candidates["weight_source"] = candidates["flow_weight_source"]
    candidates.to_csv(f"{OUT}/lr_overlap_c4_candidates.csv", index=False)
    selected = candidates[candidates["selected"]].copy()
    selected.to_csv(f"{OUT}/lr_overlap_c4_selected_targets.csv", index=False)
    coverage = (
        candidates.groupby("program", as_index=False)
        .agg(mapped_weight=("cci_weight", "sum"))
        .merge(
            selected.groupby("program", as_index=False).agg(
                selected_weight=("cci_weight", "sum"),
                selected_rows=("flow_row_id", "size"),
            ),
            on="program",
            how="left",
        )
    )
    coverage["selected_fraction"] = (
        coverage["selected_weight"] / coverage["mapped_weight"]
    )
    coverage.to_csv(f"{OUT}/lr_overlap_top10_coverage.csv", index=False)
    return candidates, selected



def build_river_rows(selected: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in selected.iterrows():
        genes = [g for g in str(row["overlap_genes"]).split(";") if g]
        if not genes:
            genes = ["no-overlap"]
        for gene in genes:
            rec = row.to_dict()
            rec["overlap_gene"] = gene
            rec["weight"] = float(row["cci_weight"]) / len(genes)
            rows.append(rec)
    river = pd.DataFrame(rows)
    river["program_label"] = river["program"].map(PROGRAM_LABELS)
    river = river.sort_values(["program_order", "rank_within_program", "overlap_gene"]).reset_index(drop=True)
    river["gi"] = np.arange(len(river))
    river.to_csv(f"{OUT}/lr_overlap_river_flows.csv", index=False)
    return river




def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    flow = load_c4_lr_flow()
    flow.to_csv(f"{OUT}/lr_overlap_active_c4_lr_input.csv", index=False)
    mapping = build_active_lr_mapping(flow)
    candidates, selected = select_c4_targets(mapping)
    river = build_river_rows(selected)

    method = {
        "input": os.path.basename(LR_AWARE_FLOW_META),
        "receiver": "Hep_C4",
        "spot_scope": "tumor sections; ST10c maxCM == CM02",
        "lr_scope": "active COMMOT LR rows with positive expression-gated score",
        "program_gene_source": os.path.basename(PROGRAM_GENES),
        "target_rule": "official pathway genes intersected with Hep_C4 expressed background; no DE filter",
        "overlap_rule": "active ligand/receptor subunit genes intersect program genes",
        "interpretation": "semantic LR gene-set overlap only; not downstream program activation or causality",
        "mechanism_pathways": {
            program: sorted(pathways)
            for program, pathways in MECHANISM_PATHWAYS.items()
        },
        "multi_program_rule": "equal split across matched programs; no duplicate score",
        "selection": f"top {SELECT_TOP_N} active LR-sender rows per program",
        "missing_program_rule": "within each predefined pathway whitelist, omit gene sets without active LR overlap; no imputation",
        "replicate_weighting": "B11 patient-balanced mean interaction-score share",
        "nichenet_used": False,
    }
    with open(
        f"{OUT}/lr_overlap_program_mapping_method.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(method, handle, ensure_ascii=False, indent=2)

    print("\nLR-overlap selected targets:")
    cols = [
        "program",
        "sender_label",
        "pathway",
        "ligand",
        "receptor",
        "overlap_genes",
        "program_match_count",
        "overlap_n",
        "cci_weight",
        "support_fraction",
    ]
    print(selected[cols].sort_values(["program", "cci_weight"], ascending=[True, False]).to_string(index=False))
    print("\nTop overlap mapping per program:")
    print(
        mapping.groupby("program", group_keys=False)
        .head(8)[
            [
                "program",
                "pathway",
                "ligand",
                "receptor",
                "overlap_n",
                "overlap_genes",
                "cci_weight",
            ]
        ]
        .to_string(index=False)
    )
    print(f"\nSaved LR-overlap outputs to {OUT}")


if __name__ == "__main__":
    main()
