"""LR-expression-aware decomposition of spot-level COMMOT communication."""

from __future__ import annotations

import gc
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


@dataclass(frozen=True)
class PairwiseConfig:
    commot_dir: Path
    output_dir: Path
    sections: tuple[str, ...]
    cm02_types: tuple[str, ...]
    display_labels: dict[str, str]
    mean_profile_file: Path
    pct_profile_file: Path
    min_detection_fraction: float = 0.05
    top_pathways_per_edge: int = 5


@dataclass(frozen=True)
class ReferenceProfiles:
    mean_cpm: pd.DataFrame
    detection_fraction: pd.DataFrame


@dataclass(frozen=True)
class PairwiseResults:
    section_lr_edges: pd.DataFrame
    patient_lr_edges: pd.DataFrame
    lr_edges: pd.DataFrame
    pathway_edges: pd.DataFrame
    edges: pd.DataFrame
    pathway_summary: pd.DataFrame
    matrix: pd.DataFrame
    prior_table: pd.DataFrame


EXPECTED_CACHE_VERSION = "cm02_indexfix_heteromeric_min_v2"
EXPECTED_HETEROMERIC_RULE = "min"


def load_reference_profiles(config: PairwiseConfig) -> ReferenceProfiles:
    mean_cpm = pd.read_csv(config.mean_profile_file, index_col="gene")
    detection = pd.read_csv(config.pct_profile_file, index_col="gene")
    common_genes = mean_cpm.index.intersection(detection.index)
    common_types = mean_cpm.columns.intersection(detection.columns)
    return ReferenceProfiles(
        mean_cpm=mean_cpm.loc[common_genes, common_types].astype(float),
        detection_fraction=detection.loc[common_genes, common_types].astype(float),
    )


def split_complex(complex_name: str) -> list[str]:
    return [gene for gene in str(complex_name).split("_") if gene]


def patient_from_section(section: str) -> str:
    match = re.match(r"^(HCC-\d+)", section)
    if match is None:
        raise ValueError(f"Cannot derive patient from section: {section}")
    return match.group(1)


def complex_expression_prior(
    complex_name: str,
    profiles: ReferenceProfiles,
    min_detection_fraction: float,
) -> pd.Series:
    genes = split_complex(complex_name)
    if not genes or any(gene not in profiles.mean_cpm.index for gene in genes):
        return pd.Series(0.0, index=profiles.mean_cpm.columns)

    expression = profiles.mean_cpm.loc[genes].clip(lower=0.0)
    detection = profiles.detection_fraction.loc[genes]
    valid = (detection >= min_detection_fraction).all(axis=0)
    positive = (expression > 0).all(axis=0)

    limiting_subunit = expression.min(axis=0)
    prior = limiting_subunit.where(valid & positive, 0.0)
    max_value = float(prior.max())
    if max_value > 0:
        prior = prior / max_value
    return prior.astype(float)


def build_complex_priors(
    ligrec: pd.DataFrame,
    profiles: ReferenceProfiles,
    config: PairwiseConfig,
) -> dict[str, pd.Series]:
    complexes = sorted(
        set(ligrec["ligand"].astype(str)).union(ligrec["receptor"].astype(str))
    )
    return {
        complex_name: complex_expression_prior(
            complex_name,
            profiles,
            config.min_detection_fraction,
        )
        for complex_name in complexes
    }


def gated_weight_matrix(
    weights: np.ndarray,
    prior: np.ndarray,
    cm02_indices: list[int],
) -> np.ndarray:
    return weights[:, cm02_indices] * prior[np.newaxis, cm02_indices]


def section_lr_flow(
    section: str,
    config: PairwiseConfig,
    profiles: ReferenceProfiles,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    adata = sc.read_h5ad(config.commot_dir / f"commot_{section}.h5ad")
    cache_version = str(adata.uns.get("cci_cache_version", ""))
    heteromeric = bool(adata.uns.get("commot_heteromeric", False))
    heteromeric_rule = str(adata.uns.get("commot_heteromeric_rule", ""))
    if cache_version != EXPECTED_CACHE_VERSION:
        raise RuntimeError(
            f"{section}: expected cache {EXPECTED_CACHE_VERSION}, got {cache_version!r}"
        )
    if not heteromeric or heteromeric_rule != EXPECTED_HETEROMERIC_RULE:
        raise RuntimeError(
            f"{section}: expected heteromeric=True/rule=min, got "
            f"{heteromeric}/{heteromeric_rule!r}"
        )
    if "rctd_all" not in adata.obsm:
        raise KeyError(f"{section}: corrected cache lacks rctd_all")

    weight_columns = list(adata.uns["rctd_all_cols"])
    weights = pd.DataFrame(
        adata.obsm["rctd_all"],
        index=adata.obs_names,
        columns=weight_columns,
    ).astype(float)
    missing_cm02 = set(config.cm02_types).difference(weight_columns)
    if missing_cm02:
        raise KeyError(f"{section}: missing CM02 RCTD types {sorted(missing_cm02)}")

    profile_types = set(profiles.mean_cpm.columns)
    unsupported_types = [celltype for celltype in weight_columns if celltype not in profile_types]
    if unsupported_types:
        weights = weights.drop(columns=unsupported_types)
        weight_columns = list(weights.columns)

    cm02_indices = [weight_columns.index(celltype) for celltype in config.cm02_types]
    weight_array = weights.to_numpy(dtype=float)
    ligrec = pd.DataFrame(adata.uns["commot-CellChat-info"]["df_ligrec"])
    ligrec = ligrec[["ligand", "receptor", "pathway"]].drop_duplicates()
    ligrec = ligrec.astype(str)
    priors = build_complex_priors(ligrec, profiles, config)

    prior_rows = []
    flow_rows = []
    for lr_row in ligrec.itertuples(index=False):
        matrix_key = f"commot-CellChat-{lr_row.ligand}-{lr_row.receptor}"
        if matrix_key not in adata.obsp:
            continue

        ligand_prior = priors[lr_row.ligand].reindex(weight_columns).fillna(0.0)
        receptor_prior = priors[lr_row.receptor].reindex(weight_columns).fillna(0.0)
        if ligand_prior.loc[list(config.cm02_types)].max() <= 0:
            continue
        if receptor_prior.loc[list(config.cm02_types)].max() <= 0:
            continue

        sender_projection = gated_weight_matrix(
            weight_array,
            ligand_prior.to_numpy(dtype=float),
            cm02_indices,
        )
        receiver_projection = gated_weight_matrix(
            weight_array,
            receptor_prior.to_numpy(dtype=float),
            cm02_indices,
        )
        communication = adata.obsp[matrix_key]
        flow = sender_projection.T.dot(communication.dot(receiver_projection))

        for celltype in weight_columns:
            prior_rows.append(
                {
                    "section": section,
                    "ligand": lr_row.ligand,
                    "receptor": lr_row.receptor,
                    "pathway": lr_row.pathway,
                    "celltype": celltype,
                    "ligand_prior": float(ligand_prior.loc[celltype]),
                    "receptor_prior": float(receptor_prior.loc[celltype]),
                }
            )

        for sender_idx, sender in enumerate(config.cm02_types):
            for receiver_idx, receiver in enumerate(config.cm02_types):
                if sender == receiver:
                    continue
                value = float(flow[sender_idx, receiver_idx])
                if value <= 0:
                    continue
                flow_rows.append(
                    {
                        "section": section,
                        "ligand": lr_row.ligand,
                        "receptor": lr_row.receptor,
                        "pathway": lr_row.pathway,
                        "sender": sender,
                        "sender_label": config.display_labels[sender],
                        "receiver": receiver,
                        "receiver_label": config.display_labels[receiver],
                        "projected_flow": value,
                    }
                )

    del adata
    gc.collect()
    edges = pd.DataFrame(flow_rows)
    if edges.empty:
        raise RuntimeError(f"{section}: no LR-aware pairwise flows")
    section_total = float(edges["projected_flow"].sum())
    edges["section_flow_share"] = edges["projected_flow"] / section_total
    priors_out = pd.DataFrame(prior_rows).drop_duplicates()
    return edges, priors_out


def aggregate_pairwise_results(
    section_edges: pd.DataFrame,
    prior_table: pd.DataFrame,
    config: PairwiseConfig,
) -> PairwiseResults:
    section_edges = section_edges.copy()
    section_edges["patient"] = section_edges["section"].map(patient_from_section)
    section_edges["lr_pair"] = (
        section_edges["ligand"].astype(str) + "->" + section_edges["receptor"].astype(str)
    )
    patients = list(dict.fromkeys(patient_from_section(section) for section in config.sections))
    section_counts = pd.Series(
        [patient_from_section(section) for section in config.sections]
    ).value_counts()
    lr_keys = [
        "ligand",
        "receptor",
        "pathway",
        "sender",
        "sender_label",
        "receiver",
        "receiver_label",
    ]
    patient_lr_edges = (
        section_edges.groupby(["patient", *lr_keys], as_index=False)
        .agg(
            section_flow_share_sum=("section_flow_share", "sum"),
            projected_flow_sum=("projected_flow", "sum"),
            support_sections=("section", "nunique"),
        )
    )
    patient_lr_edges["patient_section_count"] = patient_lr_edges["patient"].map(
        section_counts
    )
    patient_lr_edges["patient_flow_share"] = (
        patient_lr_edges["section_flow_share_sum"]
        / patient_lr_edges["patient_section_count"]
    )
    patient_lr_edges["lr_pair"] = (
        patient_lr_edges["ligand"].astype(str)
        + "->"
        + patient_lr_edges["receptor"].astype(str)
    )

    lr_edges = (
        patient_lr_edges.groupby(lr_keys, as_index=False)
        .agg(
            patient_flow_share_sum=("patient_flow_share", "sum"),
            projected_flow_sum=("projected_flow_sum", "sum"),
            support_sections=("support_sections", "sum"),
            support_patients=("patient", "nunique"),
        )
    )
    lr_edges["mean_patient_flow_share"] = (
        lr_edges["patient_flow_share_sum"] / len(patients)
    )
    lr_edges["support_fraction"] = lr_edges["support_sections"] / len(config.sections)
    lr_edges["patient_support_fraction"] = (
        lr_edges["support_patients"] / len(patients)
    )
    lr_edges = lr_edges.sort_values("mean_patient_flow_share", ascending=False)

    pathway_keys = [
        "pathway",
        "sender",
        "sender_label",
        "receiver",
        "receiver_label",
    ]
    patient_pathway_edges = (
        section_edges.groupby(["patient", *pathway_keys], as_index=False)
        .agg(
            section_flow_share_sum=("section_flow_share", "sum"),
            support_sections=("section", "nunique"),
        )
    )
    patient_pathway_edges["patient_section_count"] = patient_pathway_edges[
        "patient"
    ].map(section_counts)
    patient_pathway_edges["patient_flow_share"] = (
        patient_pathway_edges["section_flow_share_sum"]
        / patient_pathway_edges["patient_section_count"]
    )
    pathway_stats = (
        section_edges.groupby(
            pathway_keys,
            as_index=False,
        )
        .agg(
            projected_flow_sum=("projected_flow", "sum"),
            n_lr_pairs=("lr_pair", "nunique"),
            support_sections=("section", "nunique"),
        )
    )
    pathway_edges = (
        patient_pathway_edges.groupby(pathway_keys, as_index=False)
        .agg(
            patient_flow_share_sum=("patient_flow_share", "sum"),
            support_patients=("patient", "nunique"),
        )
        .merge(pathway_stats, on=pathway_keys, how="left", validate="one_to_one")
    )
    pathway_edges["mean_patient_flow_share"] = (
        pathway_edges["patient_flow_share_sum"] / len(patients)
    )
    pathway_edges["support_fraction"] = (
        pathway_edges["support_sections"] / len(config.sections)
    )
    pathway_edges["patient_support_fraction"] = (
        pathway_edges["support_patients"] / len(patients)
    )
    pathway_edges = pathway_edges.sort_values("mean_patient_flow_share", ascending=False)

    edge_rows = []
    for keys, group in pathway_edges.groupby(
        ["sender", "sender_label", "receiver", "receiver_label"], sort=False
    ):
        sender, sender_label, receiver, receiver_label = keys
        ordered = group.sort_values("mean_patient_flow_share", ascending=False)
        top = ordered.head(config.top_pathways_per_edge)
        edge_rows.append(
            {
                "sender": sender,
                "sender_label": sender_label,
                "receiver": receiver,
                "receiver_label": receiver_label,
                "mean_patient_flow_share": float(group["mean_patient_flow_share"].sum()),
                "projected_flow_sum": float(group["projected_flow_sum"].sum()),
                "n_active_pathways": int((group["mean_patient_flow_share"] > 0).sum()),
                "top_pathways": ";".join(
                    f"{row.pathway}:{row.mean_patient_flow_share:.3g}"
                    for row in top.itertuples()
                ),
            }
        )
    edges = pd.DataFrame(edge_rows).sort_values("mean_patient_flow_share", ascending=False)

    pathway_summary = (
        pathway_edges.groupby("pathway", as_index=False)
        .agg(
            mean_patient_flow_share=("mean_patient_flow_share", "sum"),
            projected_flow_sum=("projected_flow_sum", "sum"),
            n_active_edges=("mean_patient_flow_share", lambda values: int((values > 0).sum())),
            n_lr_edge_rows=("n_lr_pairs", "sum"),
        )
        .sort_values("mean_patient_flow_share", ascending=False)
    )

    matrix = pd.DataFrame(0.0, index=config.cm02_types, columns=config.cm02_types)
    for row in edges.itertuples():
        matrix.loc[row.sender, row.receiver] = row.mean_patient_flow_share

    return PairwiseResults(
        section_lr_edges=section_edges,
        patient_lr_edges=patient_lr_edges,
        lr_edges=lr_edges,
        pathway_edges=pathway_edges,
        edges=edges,
        pathway_summary=pathway_summary,
        matrix=matrix,
        prior_table=prior_table,
    )


def run_pairwise_analysis(config: PairwiseConfig) -> PairwiseResults:
    profiles = load_reference_profiles(config)
    edge_frames = []
    prior_frames = []
    for section in config.sections:
        edges, priors = section_lr_flow(section, config, profiles)
        edge_frames.append(edges)
        prior_frames.append(priors)
        print(
            f"{section}: {len(edges):,} positive LR-celltype edges; "
            f"projected flow={edges['projected_flow'].sum():.4g}"
        )
    return aggregate_pairwise_results(
        pd.concat(edge_frames, ignore_index=True),
        pd.concat(prior_frames, ignore_index=True).drop_duplicates(),
        config,
    )


def save_pairwise_results(results: PairwiseResults, config: PairwiseConfig) -> None:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    results.section_lr_edges.to_csv(
        config.output_dir / "true_pairwise_lr_edges_by_section.csv", index=False
    )
    results.patient_lr_edges.to_csv(
        config.output_dir / "true_pairwise_lr_edges_by_patient.csv", index=False
    )
    results.lr_edges.to_csv(config.output_dir / "true_pairwise_lr_edges.csv", index=False)
    results.pathway_edges.to_csv(
        config.output_dir / "true_pairwise_pathway_edges.csv", index=False
    )
    results.edges.to_csv(config.output_dir / "true_pairwise_edges.csv", index=False)
    results.pathway_summary.to_csv(
        config.output_dir / "true_pairwise_pathway_summary.csv", index=False
    )
    results.matrix.to_csv(config.output_dir / "true_pairwise_matrix_flow_share.csv")
    results.prior_table.to_csv(
        config.output_dir / "true_pairwise_lr_celltype_priors.csv", index=False
    )
    section_totals = (
        results.section_lr_edges.groupby("section", as_index=False)
        .agg(
            projected_interaction_score_sum=("projected_flow", "sum"),
            positive_lr_celltype_edges=("projected_flow", "size"),
            active_lr_pairs=(
                "lr_pair",
                "nunique",
            ),
        )
        .sort_values("section")
    )
    section_totals.to_csv(
        config.output_dir / "true_pairwise_section_totals.csv", index=False
    )
    patient_totals = (
        results.patient_lr_edges.groupby("patient", as_index=False)
        .agg(
            patient_flow_share_sum=("patient_flow_share", "sum"),
            projected_interaction_score_sum=("projected_flow_sum", "sum"),
            active_lr_pairs=(
                "lr_pair",
                "nunique",
            ),
        )
        .sort_values("patient")
    )
    patient_totals.to_csv(
        config.output_dir / "true_pairwise_patient_totals.csv", index=False
    )
    method = {
        "spot_selection": "ST10c maxCM == CM02 with corrected barcode/1-based index mapping",
        "commot_cache_version": EXPECTED_CACHE_VERSION,
        "commot_heteromeric": True,
        "commot_heteromeric_rule": "minimum subunit expression",
        "celltype_assignment": "continuous RCTD weights; no dominant label",
        "lr_assignment": "expression-gated RCTD projection without spot-wise renormalization",
        "complex_rule": "minimum subunit expression; every subunit detection fraction >= threshold",
        "min_detection_fraction": config.min_detection_fraction,
        "replicate_weighting": "within-section share; sections averaged within patient; five patients equally weighted",
        "main_edge_weight": "mean patient flow share of LR-expression-gated COMMOT interaction score",
        "share_denominator": "positive non-self CM02-to-CM02 projected scores within each tumor section",
        "interpretation": "conditional CM02 communication topology, not absolute molecular flux",
        "abundance_division": False,
        "z_scale": False,
        "sections": list(config.sections),
        "patients": list(
            dict.fromkeys(patient_from_section(section) for section in config.sections)
        ),
    }
    with open(config.output_dir / "true_pairwise_method.json", "w", encoding="utf-8") as handle:
        json.dump(method, handle, indent=2, ensure_ascii=False)
