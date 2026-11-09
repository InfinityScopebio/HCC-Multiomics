#!/usr/bin/env python3
"""Contextualize current COMMOT-supported Hep_C4 receptors with MAPK CORNETO."""

from __future__ import annotations

import json
import logging
import os
import sys
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path

import corneto
import numpy as np
import pandas as pd

from corneto_mapk_pkn import build_mapk_pkn, split_genes


LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class CornetoConfig:
    root: Path
    expression_pct_threshold: float = 5.0
    measurement_expression_pct_threshold: float = 10.0
    beta_weight: float = 0.20

    @property
    def output_dir(self) -> Path:
        return self.root / "hcc_st"

    @property
    def mapk_mapping(self) -> Path:
        return self.root / "hcc_st/CCI_nichenet_commot/lr_overlap_c4_selected_targets.csv"

    @property
    def expression(self) -> Path:
        return self.root / "hcc_st/step5_HepC4_HCC_expression.csv"

    @property
    def expressed_targets(self) -> Path:
        return self.root / "hcc_st/CCI_nichenet/target_sets_4pathway_expressed_long.csv"


def select_receptors(config: CornetoConfig) -> tuple[dict[str, int], pd.DataFrame]:
    """Select PKN inputs from the exact MAPK rows displayed by the B10 riverplot."""
    selected = pd.read_csv(config.mapk_mapping)
    selected_flag = selected["selected"].astype(str).str.lower().eq("true")
    selected = selected[(selected["program"] == "MAPK") & selected_flag].copy()
    selected["receptor_anchor"] = selected.apply(
        lambda row: ";".join(
            sorted(split_genes(row["receptor"]) & split_genes(row["overlap_genes"]))
        ),
        axis=1,
    )
    if (selected["receptor_anchor"].str.count(";") > 0).any() or (
        selected["receptor_anchor"] == ""
    ).any():
        raise ValueError("MAPK river rows do not map one-to-one to receptor anchors")
    audit = (
        selected.groupby("receptor_anchor", as_index=False)
        .agg(
            assigned_cci_weight=("cci_weight", "sum"),
            n_lr_rows=("flow_row_id", "nunique"),
            n_senders=("sender", "nunique"),
            minimum_patient_support=("patient_support_fraction", "min"),
        )
        .sort_values("assigned_cci_weight", ascending=False)
    )
    supported = set(audit["receptor_anchor"].astype(str))
    ordered = [
        anchor for anchor in ("FGFR2", "TGFBR1", "TGFBR2") if anchor in supported
    ]
    if not ordered:
        raise ValueError("No MAPK riverplot receptor anchors were found")
    return {receptor: 1 for receptor in ordered}, audit


def filter_pkn_by_expression(
    config: CornetoConfig,
    pkn_edges: list[tuple[str, int, str]],
) -> tuple[list[tuple[str, int, str]], pd.DataFrame, list[str]]:
    """Require both PKN endpoints in at least 5% of HCC Hep_C4 cells."""
    expression = pd.read_csv(config.expression)
    value_column = "Hep_C4_HCC_pct_positive"
    required_columns = {"gene", value_column}
    if not required_columns.issubset(expression.columns):
        raise ValueError(f"Expression table lacks columns: {sorted(required_columns)}")
    expression_by_gene = expression.set_index("gene")[value_column]
    pkn_nodes = {source for source, _, _ in pkn_edges} | {
        target for _, _, target in pkn_edges
    }
    missing = sorted(pkn_nodes.difference(expression_by_gene.index))
    if missing:
        raise ValueError(f"HCC Hep_C4 expression is missing PKN nodes: {missing}")
    retained = [
        edge
        for edge in pkn_edges
        if expression_by_gene.loc[edge[0]] >= config.expression_pct_threshold
        and expression_by_gene.loc[edge[2]] >= config.expression_pct_threshold
    ]
    retained_nodes = {source for source, _, _ in retained} | {
        target for _, _, target in retained
    }
    dropped_nodes = sorted(pkn_nodes.difference(retained_nodes))
    return retained, expression, dropped_nodes


def select_measurements(
    config: CornetoConfig,
    expression: pd.DataFrame,
    pkn_edges: list[tuple[str, int, str]],
) -> tuple[dict[str, int], pd.DataFrame]:
    """Select current no-DE, expressed MAPK targets at the PKN terminus."""
    candidates = ("MYC", "FOS", "JUN", "ETS1", "ELK1", "CREB1")
    target_table = pd.read_csv(config.expressed_targets)
    required_columns = {"program", "gene", "gene_set_name", "target_filter"}
    if not required_columns.issubset(target_table.columns):
        raise ValueError(f"Target table lacks columns: {sorted(required_columns)}")
    mapk_targets = target_table[target_table["program"] == "MAPK"].copy()
    if mapk_targets.empty:
        raise ValueError("Current expressed-target table contains no MAPK genes")
    target_filters = set(mapk_targets["target_filter"].astype(str))
    if target_filters != {
        "official_pathway_intersect_hepc4_expressed_background_no_DE"
    }:
        raise ValueError(f"Unexpected MAPK target filters: {sorted(target_filters)}")
    current_mapk_genes = set(mapk_targets["gene"].astype(str))
    expression_indexed = expression.set_index("gene")
    rows: list[dict[str, float | str | bool]] = []
    retained_nodes = {source for source, _, _ in pkn_edges} | {
        target for _, _, target in pkn_edges
    }
    for gene in candidates:
        pct = float(expression_indexed.at[gene, "Hep_C4_HCC_pct_positive"])
        selected = bool(
            gene in retained_nodes
            and gene in current_mapk_genes
            and pct >= config.measurement_expression_pct_threshold
        )
        rows.append(
            {
                "gene": gene,
                "in_current_MAPK_expressed_target_set": gene in current_mapk_genes,
                "Hep_C4_HCC_pct_positive": pct,
                "selected": selected,
            }
        )
    audit = pd.DataFrame(rows)
    selected_genes = audit.loc[audit["selected"], "gene"].tolist()
    if not selected_genes:
        raise ValueError("No current expressed MAPK targets passed CORNETO filters")
    return {gene: 1 for gene in selected_genes}, audit


def extract_selected_edges(
    problem: object,
    graph: object,
    pkn_edges: list[tuple[str, int, str]],
    receptor: str,
) -> tuple[pd.DataFrame, set[str]]:
    """Extract PKN edges carrying nonzero CARNIVAL flow."""
    flow = np.asarray(problem.symbols["_flow"].value, dtype=float).reshape(-1)
    activities = np.asarray(
        problem.symbols["species_activated_c0"].value, dtype=float
    ).reshape(-1)
    graph_nodes = list(graph.V)
    activated = {
        graph_nodes[index]
        for index, value in enumerate(activities)
        if value > 0.5 and not str(graph_nodes[index]).startswith("_")
    }
    pkn_lookup = defaultdict(list)
    for source, sign, target in pkn_edges:
        pkn_lookup[(source, target)].append(sign)
    rows: list[dict[str, str | int | float]] = []
    for index, (source_set, target_set) in enumerate(graph.E):
        if index >= len(flow) or flow[index] <= 1e-6:
            continue
        sources = [node for node in source_set if not str(node).startswith("_")]
        targets = [node for node in target_set if not str(node).startswith("_")]
        if len(sources) != 1 or len(targets) != 1:
            continue
        source, target = str(sources[0]), str(targets[0])
        for sign in pkn_lookup.get((source, target), []):
            rows.append(
                {
                    "receptor": receptor,
                    "source": source,
                    "target": target,
                    "sign": sign,
                    "flow": float(flow[index]),
                }
            )
    columns = ["receptor", "source", "target", "sign", "flow"]
    selected = pd.DataFrame(rows, columns=columns).drop_duplicates(
        subset=["receptor", "source", "target", "sign"]
    )
    return selected.reset_index(drop=True), activated


def run_carnival(
    perturbations: dict[str, int],
    measurements: dict[str, int],
    pkn_edges: list[tuple[str, int, str]],
    beta_weight: float,
    receptor_label: str,
) -> tuple[pd.DataFrame, set[str]]:
    """Run one CARNIVAL model and return its selected PKN edges."""
    problem, graph = corneto.methods.runVanillaCarnival(
        perturbations=perturbations,
        measurements=measurements,
        priorKnowledgeNetwork=pkn_edges,
        betaWeight=beta_weight,
        solver="SCIPY",
        backend_options={
            "scipy_options": {"random_seed": 0, "threads": 1},
        },
        verbose=False,
    )
    return extract_selected_edges(problem, graph, pkn_edges, receptor_label)


def reachable_targets(
    receptor: str,
    edges: pd.DataFrame,
    measurements: dict[str, int],
) -> tuple[list[dict[str, str | int]], set[str]]:
    """Audit which measured targets are reachable in a selected subnetwork."""
    adjacency: dict[str, list[str]] = defaultdict(list)
    for row in edges.itertuples(index=False):
        if int(row.sign) > 0:
            adjacency[str(row.source)].append(str(row.target))
    seen = {receptor}
    queue: deque[str] = deque([receptor])
    while queue:
        source = queue.popleft()
        for target in adjacency.get(source, []):
            if target not in seen:
                seen.add(target)
                queue.append(target)
    rows = [
        {"receptor": receptor, "target": target, "reached": int(target in seen)}
        for target in measurements
    ]
    return rows, seen


def save_method_audit(
    config: CornetoConfig,
    perturbations: dict[str, int],
    measurements: dict[str, int],
    pkn_before: int,
    pkn_after: int,
    dropped_nodes: list[str],
) -> None:
    """Record the current CORNETO input and interpretation contract."""
    method = {
        "analysis": "COMMOT-anchored MAPK CORNETO contextualization",
        "cci_input": str(config.mapk_mapping.relative_to(config.root)),
        "cci_scope": (
            "B10 MAPK riverplot top-10 rows derived from B11 tumor sections; "
            "maxCM == CM02 spots; receiver Hep_C4"
        ),
        "receptor_selection": (
            "MAPK rows selected by active LR-gene overlap with the current official "
            "MAPK x HCC Hep_C4 expressed-background target set; receptor overlap gene "
            "used as the PKN anchor"
        ),
        "perturbations": list(perturbations),
        "expression_input": str(config.expression.relative_to(config.root)),
        "expression_scope": "947 HCC Hep_C4 cells",
        "expression_pct_threshold": config.expression_pct_threshold,
        "measurement_input": str(config.expressed_targets.relative_to(config.root)),
        "measurement_rule": (
            "curated PKN ERK-terminal candidate intersected with the current "
            "KEGG MAPK pathway x HCC Hep_C4 expressed-background target set; "
            f"no DE; HCC Hep_C4 positive >= "
            f"{config.measurement_expression_pct_threshold}%"
        ),
        "measurements": list(measurements),
        "pkn_edges_before_expression_filter": pkn_before,
        "pkn_edges_after_expression_filter": pkn_after,
        "dropped_pkn_nodes": dropped_nodes,
        "beta_weight": config.beta_weight,
        "python_hash_seed": 0,
        "solver_note": (
            "SCIPY/HiGHS MILP with random_seed=0 and threads=1; "
            "PYTHONHASHSEED=0 also enforces deterministic graph ordering"
        ),
        "interpretation": (
            "CORNETO selects sign-consistent routes within a curated PKN. It does not "
            "estimate CCI edge weights, prove signal direction in tissue, or establish causality."
        ),
    }
    output = config.output_dir / "step5_corneto_method.json"
    output.write_text(json.dumps(method, indent=2), encoding="utf-8")


def ensure_deterministic_hash_seed() -> None:
    """Re-exec once because CORNETO graph construction traverses Python sets."""
    if os.environ.get("PYTHONHASHSEED") == "0":
        return
    environment = dict(os.environ)
    environment["PYTHONHASHSEED"] = "0"
    os.execve(sys.executable, [sys.executable, *sys.argv], environment)


def main() -> None:
    ensure_deterministic_hash_seed()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config = CornetoConfig(root=Path(__file__).resolve().parents[1])
    raw_pkn = build_mapk_pkn()
    perturbations, receptor_audit = select_receptors(config)
    pkn_edges, expression, dropped_nodes = filter_pkn_by_expression(config, raw_pkn)
    measurements, measurement_audit = select_measurements(config, expression, pkn_edges)

    receptor_audit.to_csv(
        config.output_dir / "step5_corneto_input_receptors.csv", index=False
    )
    measurement_audit.to_csv(
        config.output_dir / "step5_corneto_measurements.csv", index=False
    )
    pd.DataFrame(pkn_edges, columns=["source", "sign", "target"]).to_csv(
        config.output_dir / "step5_corneto_filtered_pkn.csv", index=False
    )

    joint_edges, joint_activated = run_carnival(
        perturbations, measurements, pkn_edges, config.beta_weight, "joint"
    )
    joint_edges.to_csv(config.output_dir / "step5_corneto_active_edges.csv", index=False)
    pd.DataFrame({"node": sorted(joint_activated)}).to_csv(
        config.output_dir / "step5_corneto_activated_nodes.csv", index=False
    )

    per_receptor: list[pd.DataFrame] = []
    reachability_rows: list[dict[str, str | int]] = []
    seen_by_receptor: dict[str, set[str]] = {}
    for receptor in perturbations:
        selected, _ = run_carnival(
            {receptor: 1}, measurements, pkn_edges, config.beta_weight, receptor
        )
        reached, seen = reachable_targets(receptor, selected, measurements)
        if not all(row["reached"] for row in reached):
            missing = [row["target"] for row in reached if not row["reached"]]
            raise ValueError(f"{receptor} does not reach measurements: {missing}")
        per_receptor.append(selected)
        reachability_rows.extend(reached)
        seen_by_receptor[receptor] = seen
    per_receptor_edges = pd.concat(per_receptor, ignore_index=True)
    per_receptor_edges.to_csv(
        config.output_dir / "step5_corneto_per_receptor.csv", index=False
    )
    pd.DataFrame(reachability_rows).to_csv(
        config.output_dir / "step5_corneto_reachability.csv", index=False
    )

    measurement_nodes = set(measurements)
    convergence = defaultdict(int)
    for receptor, seen in seen_by_receptor.items():
        for node in seen.difference({receptor}).difference(measurement_nodes):
            convergence[node] += 1
    convergence_table = pd.DataFrame(
        [
            {"node": node, "n_receptors_reaching": count}
            for node, count in convergence.items()
        ]
    ).sort_values(["n_receptors_reaching", "node"], ascending=[False, True])
    convergence_table.to_csv(
        config.output_dir / "step5_corneto_convergence.csv", index=False
    )
    save_method_audit(
        config,
        perturbations,
        measurements,
        len(raw_pkn),
        len(pkn_edges),
        dropped_nodes,
    )
    LOGGER.info(
        "CORNETO complete: %d receptors, %d measurements, %d filtered PKN edges",
        len(perturbations),
        len(measurements),
        len(pkn_edges),
    )
    LOGGER.info("Measurements: %s", ", ".join(measurements))


if __name__ == "__main__":
    main()
