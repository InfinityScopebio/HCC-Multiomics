#!/usr/bin/env python3
"""Run current CORNETO twice without frozen-route fallback."""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import json
import os
import random
import shutil
import subprocess
import sys
import warnings
from pathlib import Path
from types import ModuleType
from typing import Any

import corneto
import numpy as np
import pandas as pd
from corneto import Graph
from corneto.backend._base import ProblemDef
from corneto.methods.signflow import create_flow_graph, signflow

from config import SpatialConfig, load_config


REPEAT_RUN_NAMES = ("run1", "run2")
PRIMARY_OBJECTIVE_TOLERANCE = 1e-7


def strict_boolean_series(values: pd.Series, field: str) -> pd.Series:
    """Parse persisted booleans without treating non-empty strings as true."""

    normalized = values.astype(str).str.strip().str.lower()
    invalid = ~normalized.isin({"true", "false"})
    if invalid.any():
        unexpected = sorted(normalized.loc[invalid].unique().tolist())
        raise RuntimeError(f"{field} contains invalid booleans: {unexpected}")
    return normalized.eq("true")


class _StaticSymbol:
    def __init__(self, values: list[float]) -> None:
        self.value = np.asarray(values, dtype=float)


class _StaticProblem:
    def __init__(
        self,
        activation: float | list[float],
        inhibition: float | list[float],
        flow: float | list[float],
    ) -> None:
        def vector(value: float | list[float]) -> list[float]:
            return value if isinstance(value, list) else [value]

        edge_count = len(vector(activation))
        self.symbols = {
            "reaction_sends_activation_c0": _StaticSymbol(vector(activation)),
            "reaction_sends_inhibition_c0": _StaticSymbol(vector(inhibition)),
            "species_activated_c0": _StaticSymbol([1.0, 1.0]),
            "species_inhibited_c0": _StaticSymbol([0.0, 0.0]),
            "_flow": _StaticSymbol(vector(flow)),
        }
        if not all(
            len(vector(value)) == edge_count
            for value in (inhibition, flow)
        ):
            raise ValueError("Static edge vectors must have equal lengths.")


class _StaticGraph:
    V = ("A", "B")
    E = ((frozenset({"A"}), frozenset({"B"})),)
    edge_properties = ({"interaction": 1},)


class _ParallelStaticGraph:
    V = ("A", "B")
    E = (
        (frozenset({"A"}), frozenset({"B"})),
        (frozenset({"A"}), frozenset({"B"})),
    )
    edge_properties = (
        {"interaction": 1},
        {"interaction": -1},
    )


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def copy_input(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def extract_signal_edges(
    problem: object,
    graph: object,
    pkn_edges: list[tuple[str, int, str]],
    receptor: str,
) -> tuple[pd.DataFrame, set[str]]:
    """Extract CARNIVAL signal topology without using unconstrained flow."""

    activation = np.asarray(
        problem.symbols["reaction_sends_activation_c0"].value,
        dtype=float,
    ).reshape(-1)
    inhibition = np.asarray(
        problem.symbols["reaction_sends_inhibition_c0"].value,
        dtype=float,
    ).reshape(-1)
    if np.any((activation > 0.5) & (inhibition > 0.5)):
        raise RuntimeError("A reaction carries activation and inhibition simultaneously.")

    species_activation = np.asarray(
        problem.symbols["species_activated_c0"].value,
        dtype=float,
    ).reshape(-1)
    species_inhibition = np.asarray(
        problem.symbols["species_inhibited_c0"].value,
        dtype=float,
    ).reshape(-1)
    graph_nodes = list(graph.V)
    activated = {
        str(graph_nodes[index])
        for index in range(len(graph_nodes))
        if (
            species_activation[index] > 0.5
            or species_inhibition[index] > 0.5
        )
        and not str(graph_nodes[index]).startswith("_")
    }

    pkn_set = {
        (str(source), int(sign), str(target))
        for source, sign, target in pkn_edges
    }

    rows: list[dict[str, str | int]] = []
    for index, (source_set, target_set) in enumerate(graph.E):
        sends_activation = activation[index] > 0.5
        sends_inhibition = inhibition[index] > 0.5
        if not sends_activation and not sends_inhibition:
            continue
        sources = [node for node in source_set if not str(node).startswith("_")]
        targets = [node for node in target_set if not str(node).startswith("_")]
        if len(sources) != 1 or len(targets) != 1:
            continue
        source, target = str(sources[0]), str(targets[0])
        interaction = int(
            graph.edge_properties[index].get("interaction", 0)
        )
        if (source, interaction, target) not in pkn_set:
            continue
        signal_state = 1 if sends_activation else -1
        rows.append(
            {
                "receptor": receptor,
                "source": source,
                "target": target,
                "sign": interaction,
                "signal_state": signal_state,
            }
        )
    columns = ["receptor", "source", "target", "sign", "signal_state"]
    selected = pd.DataFrame(rows, columns=columns).drop_duplicates()
    selected = selected.sort_values(columns, kind="mergesort").reset_index(drop=True)
    return selected, activated


def signal_flow_diagnostic(
    problem: object,
    graph: object,
    pkn_edges: list[tuple[str, int, str]],
    receptor: str,
) -> dict[str, str | int]:
    """Quantify why positive flow is not equivalent to selected signal."""

    flow = np.asarray(problem.symbols["_flow"].value, dtype=float).reshape(-1)
    activation = np.asarray(
        problem.symbols["reaction_sends_activation_c0"].value,
        dtype=float,
    ).reshape(-1)
    inhibition = np.asarray(
        problem.symbols["reaction_sends_inhibition_c0"].value,
        dtype=float,
    ).reshape(-1)
    pkn_set = {
        (str(source), int(sign), str(target))
        for source, sign, target in pkn_edges
    }

    flow_keys: set[tuple[str, str, int]] = set()
    signal_keys: set[tuple[str, str, int]] = set()
    for index, (source_set, target_set) in enumerate(graph.E):
        sources = [node for node in source_set if not str(node).startswith("_")]
        targets = [node for node in target_set if not str(node).startswith("_")]
        if len(sources) != 1 or len(targets) != 1:
            continue
        source, target = str(sources[0]), str(targets[0])
        interaction = int(
            graph.edge_properties[index].get("interaction", 0)
        )
        key = (source, target, interaction)
        if (source, interaction, target) not in pkn_set:
            continue
        if flow[index] > 1e-6:
            flow_keys.add(key)
        if activation[index] > 0.5 or inhibition[index] > 0.5:
            signal_keys.add(key)

    def serialize(keys: set[tuple[str, str, int]]) -> str:
        return ";".join(
            f"{source}>{target}[{sign:+d}]"
            for source, target, sign in sorted(keys)
        )

    return {
        "receptor": receptor,
        "positive_flow_edges": len(flow_keys),
        "signal_edges": len(signal_keys),
        "flow_only_edges": len(flow_keys - signal_keys),
        "signal_without_flow_edges": len(signal_keys - flow_keys),
        "flow_only_edge_keys": serialize(flow_keys - signal_keys),
        "signal_without_flow_edge_keys": serialize(signal_keys - flow_keys),
    }


def validate_signal_extractor(config: SpatialConfig) -> None:
    """Regression-test that route extraction is signal-only, never flow-based."""

    cases = (
        ("positive_flow_without_signal", 0.0, 0.0, 1.0, 0, None),
        ("activation_with_zero_flow", 1.0, 0.0, 0.0, 1, 1),
        ("inhibition_with_zero_flow", 0.0, 1.0, 0.0, 1, -1),
    )
    rows: list[dict[str, str | int | float | None]] = []
    for name, activation, inhibition, flow, expected_edges, expected_state in cases:
        selected, _ = extract_signal_edges(
            _StaticProblem(activation, inhibition, flow),
            _StaticGraph(),
            [("A", 1, "B")],
            "R",
        )
        observed_state = (
            int(selected.iloc[0]["signal_state"])
            if len(selected) == 1
            else None
        )
        passed = len(selected) == expected_edges and observed_state == expected_state
        rows.append(
            {
                "case": name,
                "activation_indicator": activation,
                "inhibition_indicator": inhibition,
                "flow": flow,
                "expected_edges": expected_edges,
                "observed_edges": len(selected),
                "expected_signal_state": expected_state,
                "observed_signal_state": observed_state,
                "status": "PASS" if passed else "FAIL",
            }
        )
    parallel_selected, _ = extract_signal_edges(
        _StaticProblem([1.0, 0.0], [0.0, 0.0], [1.0, 1.0]),
        _ParallelStaticGraph(),
        [("A", 1, "B"), ("A", -1, "B")],
        "R",
    )
    parallel_passed = (
        len(parallel_selected) == 1
        and int(parallel_selected.iloc[0]["sign"]) == 1
        and int(parallel_selected.iloc[0]["signal_state"]) == 1
    )
    rows.append(
        {
            "case": "parallel_opposite_sign_exact_triple",
            "activation_indicator": "1|0",
            "inhibition_indicator": "0|0",
            "flow": "1|1",
            "expected_edges": 1,
            "observed_edges": len(parallel_selected),
            "expected_signal_state": 1,
            "observed_signal_state": (
                int(parallel_selected.iloc[0]["signal_state"])
                if len(parallel_selected) == 1
                else None
            ),
            "status": "PASS" if parallel_passed else "FAIL",
        }
    )
    audit = pd.DataFrame(rows)
    output = config.audit_root / "spatial_corneto_signal_extractor_regression.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output, sep="\t", index=False)
    if not audit["status"].eq("PASS").all():
        raise RuntimeError("CORNETO signal-extractor regression failed.")


def solve_problem(
    problem: ProblemDef,
    config: SpatialConfig,
) -> dict[str, str | float]:
    """Solve one MILP with the audited HiGHS controls."""

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Unrecognized options detected.*",
            category=RuntimeWarning,
        )
        solved = problem.solve(
            solver="SCIPY",
            scipy_options={
                "random_seed": config.corneto_solver_seed,
                "threads": config.threads,
            },
            verbosity=0,
        )
    status = str(getattr(solved, "status", "")).lower()
    if status != "optimal":
        raise RuntimeError(f"CORNETO MILP did not solve to optimality: {status!r}.")
    solver_stats = getattr(solved, "solver_stats", None)
    extra_stats = getattr(solver_stats, "extra_stats", None)
    if not isinstance(extra_stats, dict) or "mip_gap" not in extra_stats:
        raise RuntimeError("CORNETO MILP solver did not report a MIP gap.")
    mip_gap = float(extra_stats["mip_gap"])
    if not np.isfinite(mip_gap) or not 0.0 <= mip_gap <= 1e-9:
        raise RuntimeError(f"CORNETO MILP has invalid MIP gap: {mip_gap}.")
    return {"solver_status": status, "mip_gap": mip_gap}


def build_carnival_problem(
    perturbations: dict[str, int],
    measurements: dict[str, int],
    pkn_edges: list[tuple[str, int, str]],
    beta_weight: float,
) -> tuple[ProblemDef, Graph]:
    """Build the same primary CARNIVAL model as runVanillaCarnival."""

    data: dict[str, tuple[str, int]] = {
        gene: ("P", value) for gene, value in perturbations.items()
    }
    data.update({gene: ("M", value) for gene, value in measurements.items()})
    conditions = {"c0": data}
    graph = Graph.from_sif_tuples(pkn_edges)
    flow_graph = create_flow_graph(graph, conditions)
    problem = signflow(
        flow_graph,
        conditions,
        l0_penalty_vertices=beta_weight,
    )
    return problem, flow_graph


def weighted_objective(problem: ProblemDef) -> object:
    """Return the original weighted CARNIVAL primary objective expression."""

    if not problem.objectives:
        raise RuntimeError("CARNIVAL problem has no primary objective.")
    expression = problem.objectives[0] * problem.weights[0]
    for objective, weight in zip(problem.objectives[1:], problem.weights[1:]):
        expression = expression + objective * weight
    return expression


def scalar_value(expression: object) -> float:
    """Read a solved scalar proxy expression."""

    value = np.asarray(expression.value, dtype=float)
    if value.size == 0 or not np.isfinite(value).all():
        raise RuntimeError("Solver expression has no finite value.")
    return float(value.sum())


def binary_value(expression: object) -> int:
    """Read one solved binary expression with an integrality check."""

    value = scalar_value(expression)
    rounded = int(round(value))
    if rounded not in (0, 1) or abs(value - rounded) > 1e-5:
        raise RuntimeError(f"Expected a binary solver value, observed {value}.")
    return rounded


def integer_value(
    expression: object,
    *,
    minimum: int,
    maximum: int,
) -> int:
    """Read one solved integer expression without silently rounding drift."""

    value = scalar_value(expression)
    rounded = int(round(value))
    if (
        rounded < minimum
        or rounded > maximum
        or abs(value - rounded) > 1e-5
    ):
        raise RuntimeError(
            f"Expected an integer solver value in [{minimum}, {maximum}], "
            f"observed {value}."
        )
    return rounded


def objective_problem(
    base: ProblemDef,
    extra_constraints: list[object],
    objective: object,
) -> ProblemDef:
    """Reuse one CARNIVAL model with a deterministic secondary objective."""

    return ProblemDef(
        backend=base._backend,
        constraints=[*base.constraints, *extra_constraints],
        objectives=[objective],
        weights=[1.0],
        direction=base.direction,
        graph=base._graph,
    )


def canonical_pkn_edge_indices(
    graph: Graph,
    pkn_edges: list[tuple[str, int, str]],
) -> list[tuple[tuple[str, str, int], int]]:
    """Map real PKN reactions to a stable canonical edge order."""

    pkn_set = {(source, int(sign), target) for source, sign, target in pkn_edges}
    ordered: list[tuple[tuple[str, str, int], int]] = []
    for index, (source_set, target_set) in enumerate(graph.edges):
        sources = [str(node) for node in source_set if not str(node).startswith("_")]
        targets = [str(node) for node in target_set if not str(node).startswith("_")]
        if len(sources) != 1 or len(targets) != 1:
            continue
        sign = int(graph.edge_properties[index].get("interaction", 0))
        triple = (sources[0], sign, targets[0])
        if triple in pkn_set:
            ordered.append(((sources[0], targets[0], sign), index))
    ordered.sort(key=lambda item: item[0])
    if len(ordered) != len(pkn_set):
        raise RuntimeError(
            "Canonical PKN edge mapping is not one-to-one: "
            f"{len(ordered)} graph edges for {len(pkn_set)} PKN edges."
        )
    return ordered


def run_lexicographic_carnival(
    config: SpatialConfig,
    perturbations: dict[str, int],
    measurements: dict[str, int],
    pkn_edges: list[tuple[str, int, str]],
    beta_weight: float,
    receptor_label: str,
) -> tuple[ProblemDef, Graph, dict[str, str | int | float]]:
    """Lock the primary optimum, then select one parsimonious canonical topology."""

    base, graph = build_carnival_problem(
        perturbations,
        measurements,
        pkn_edges,
        beta_weight,
    )
    solver_audits = [solve_problem(base, config)]
    primary = weighted_objective(base)
    primary_optimum = scalar_value(primary)
    constraints: list[object] = [
        primary <= primary_optimum + PRIMARY_OBJECTIVE_TOLERANCE,
        primary >= primary_optimum - PRIMARY_OBJECTIVE_TOLERANCE,
    ]

    activation = base.symbols["reaction_sends_activation_c0"]
    inhibition = base.symbols["reaction_sends_inhibition_c0"]
    ordered = canonical_pkn_edge_indices(graph, pkn_edges)
    indices = [index for _, index in ordered]
    signal_count = np.ones(len(indices)) @ (
        activation[indices] + inhibition[indices]
    )
    cardinality_problem = objective_problem(base, constraints, signal_count)
    solver_audits.append(solve_problem(cardinality_problem, config))
    minimum_signal_edges = integer_value(
        signal_count,
        minimum=0,
        maximum=len(indices),
    )
    constraints.append(signal_count == minimum_signal_edges)

    selected: dict[int, int] = {}
    for _, index in ordered:
        indicator = activation[index] + inhibition[index]
        lex_problem = objective_problem(base, constraints, -indicator)
        solver_audits.append(solve_problem(lex_problem, config))
        selected[index] = binary_value(indicator)
        constraints.append(indicator == selected[index])

    activation_states: dict[int, int] = {}
    for _, index in ordered:
        if selected[index] == 0:
            continue
        direction_problem = objective_problem(base, constraints, -activation[index])
        solver_audits.append(solve_problem(direction_problem, config))
        activation_states[index] = binary_value(activation[index])
        constraints.append(activation[index] == activation_states[index])

    final_problem = objective_problem(base, constraints, signal_count * 0.0)
    solver_audits.append(solve_problem(final_problem, config))
    final_primary = scalar_value(primary)
    final_count = integer_value(
        signal_count,
        minimum=0,
        maximum=len(indices),
    )
    if abs(final_primary - primary_optimum) > PRIMARY_OBJECTIVE_TOLERANCE:
        raise RuntimeError("Lexicographic solve changed the CARNIVAL primary optimum.")
    if final_count != minimum_signal_edges:
        raise RuntimeError("Lexicographic solve changed the secondary edge optimum.")

    audit: dict[str, str | int | float] = {
        "receptor": receptor_label,
        "primary_objective": primary_optimum,
        "final_primary_objective": final_primary,
        "primary_tolerance": PRIMARY_OBJECTIVE_TOLERANCE,
        "minimum_signal_edges": minimum_signal_edges,
        "canonical_pkn_edges": len(ordered),
        "activation_signal_edges": sum(activation_states.values()),
        "inhibition_signal_edges": len(activation_states)
        - sum(activation_states.values()),
        "solver_calls": len(solver_audits),
        "solver_statuses": "|".join(
            sorted({str(row["solver_status"]) for row in solver_audits})
        ),
        "maximum_mip_gap": max(
            float(row["mip_gap"]) for row in solver_audits
        ),
        "solver_all_optimal": all(
            row["solver_status"] == "optimal"
            and 0.0 <= float(row["mip_gap"]) <= 1e-9
            for row in solver_audits
        ),
        "master_seed": config.master_seed,
        "solver_seed": config.corneto_solver_seed,
        "threads": config.threads,
        "secondary_rule": (
            "minimize_signal_edge_count_then_canonical_first_edges_then_"
            "activation_before_inhibition"
        ),
        "status": "locked_primary_and_lexicographic_optimum",
    }
    return final_problem, graph, audit


def install_compatibility(stage: ModuleType, config: SpatialConfig) -> None:
    diagnostics: list[dict[str, str | int]] = []
    lexicographic_audits: list[dict[str, str | int | float]] = []

    def run_carnival(
        perturbations: dict[str, int],
        measurements: dict[str, int],
        pkn_edges: list[tuple[str, int, str]],
        beta_weight: float,
        receptor_label: str,
    ) -> tuple[pd.DataFrame, set[str]]:
        problem, graph, lexicographic_audit = run_lexicographic_carnival(
            config,
            perturbations,
            measurements,
            pkn_edges,
            beta_weight,
            receptor_label,
        )
        lexicographic_audits.append(lexicographic_audit)
        graph_type = type(graph)
        if not hasattr(graph_type, "V"):
            graph_type.V = property(lambda instance: instance.vertices)
        if not hasattr(graph_type, "E"):
            graph_type.E = property(lambda instance: instance.edges)
        diagnostics.append(
            signal_flow_diagnostic(
                problem,
                graph,
                pkn_edges,
                receptor_label,
            )
        )
        return extract_signal_edges(
            problem,
            graph,
            pkn_edges,
            receptor_label,
        )

    stage.run_carnival = run_carnival
    stage.ensure_deterministic_hash_seed = lambda: None
    stage.signal_flow_diagnostics = diagnostics
    stage.lexicographic_audits = lexicographic_audits


def prepare_workspace(config: SpatialConfig, run_name: str) -> Path:
    workspace = config.corneto_root / "repeat_runs" / run_name
    if workspace.exists():
        shutil.rmtree(workspace)
    hcc = workspace / "hcc_st"
    copy_input(
        config.cci_root / "lr_overlap_c4_selected_targets.csv",
        hcc / "CCI_nichenet_commot/lr_overlap_c4_selected_targets.csv",
    )
    copy_input(
        config.cci_root / "step5_HepC4_HCC_expression.csv",
        hcc / "step5_HepC4_HCC_expression.csv",
    )
    copy_input(
        config.cci_root / "target_sets_4pathway_expressed_long.csv",
        hcc / "CCI_nichenet/target_sets_4pathway_expressed_long.csv",
    )
    return workspace


def run_once(config: SpatialConfig, run_name: str) -> Path:
    random.seed(config.master_seed)
    np.random.seed(config.master_seed)
    workspace = prepare_workspace(config, run_name)
    pkn_path = config.assert_legacy_code(
        config.legacy_code_root / "corneto_mapk_pkn.py"
    )
    stage_path = config.assert_legacy_code(
        config.legacy_code_root / "run_ST_step5_corneto_pathway.py"
    )
    load_module("corneto_mapk_pkn", pkn_path)
    stage = load_module(f"fixed_corneto_{run_name}", stage_path)
    stage.__file__ = str(workspace / "hcc_st/run_ST_step5_corneto_pathway.py")
    install_compatibility(stage, config)
    stage.main()
    pd.DataFrame(stage.signal_flow_diagnostics).to_csv(
        workspace / "hcc_st/step5_corneto_signal_vs_flow_diagnostic.csv",
        index=False,
    )
    pd.DataFrame(stage.lexicographic_audits).to_csv(
        workspace / "hcc_st/step5_corneto_lexicographic_audit.csv",
        index=False,
    )
    method_path = workspace / "hcc_st/step5_corneto_method.json"
    method = json.loads(method_path.read_text(encoding="utf-8"))
    method.update(
        {
            "route_extraction": (
                "reaction_sends_activation_c0 and "
                "reaction_sends_inhibition_c0"
            ),
            "flow_used_for_topology": False,
            "flow_exclusion_reason": (
                "signal_implies_flow constrains signal <= flow indicator; "
                "positive _flow can split across equivalent routes and is not "
                "a valid selected-signal criterion"
            ),
            "repeat_process_isolation": "fresh process per complete CARNIVAL run",
            "secondary_tie_break": (
                "lock original CARNIVAL primary optimum within 1e-7; minimize "
                "real-PKN signal-edge cardinality; then canonical-first "
                "lexicographic edge selection and activation-before-inhibition"
            ),
            "python_hash_seed": config.master_seed,
            "solver_note": (
                "SCIPY/HiGHS MILP with random_seed="
                f"{config.corneto_solver_seed}, threads={config.threads}, and "
                f"PYTHONHASHSEED={config.master_seed}; each complete run uses a "
                "fresh process"
            ),
            "solver_acceptance": (
                "every MILP call must return cvxpy status=optimal and "
                "HiGHS mip_gap<=1e-9"
            ),
        }
    )
    method_path.write_text(json.dumps(method, indent=2), encoding="utf-8")
    return workspace / "hcc_st/step5_corneto_per_receptor.csv"


def route_keys(table: pd.DataFrame) -> set[tuple[str, str, str, str, str]]:
    columns = ["receptor", "source", "target", "sign", "signal_state"]
    return set(
        table[columns].astype(str).itertuples(index=False, name=None)
    )


def version_info(config: SpatialConfig) -> dict[str, Any]:
    versions: dict[str, Any] = {
        "corneto": getattr(corneto, "__version__", "unknown"),
        "python": sys.version.split()[0],
        "pythonhashseed": os.environ.get("PYTHONHASHSEED"),
        "master_seed": config.master_seed,
        "solver_seed": config.corneto_solver_seed,
        "threads": config.threads,
    }
    for module_name in ("scipy", "cvxpy", "highspy"):
        try:
            module = importlib.import_module(module_name)
            versions[module_name] = getattr(module, "__version__", "unknown")
        except ImportError:
            versions[module_name] = "unavailable"
    return versions


def configure_runtime(config: SpatialConfig) -> dict[str, str]:
    """Fix process-level RNG and thread controls for CORNETO/HiGHS."""

    environment = dict(os.environ)
    controls = {
        "PYTHONHASHSEED": str(config.master_seed),
        "OMP_NUM_THREADS": str(config.threads),
        "OPENBLAS_NUM_THREADS": str(config.threads),
        "MKL_NUM_THREADS": str(config.threads),
        "BLIS_NUM_THREADS": str(config.threads),
        "NUMEXPR_NUM_THREADS": str(config.threads),
        "RCPP_PARALLEL_NUM_THREADS": str(config.threads),
    }
    environment.update(controls)
    os.environ.update(controls)
    random.seed(config.master_seed)
    np.random.seed(config.master_seed)
    return environment


def run_isolated_worker(
    config: SpatialConfig,
    run_name: str,
    environment: dict[str, str],
) -> Path:
    """Run a complete CORNETO solve in a fresh, fixed process environment."""

    subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "--worker", run_name],
        check=True,
        cwd=config.project_root,
        env=environment,
    )
    return (
        config.corneto_root
        / "repeat_runs"
        / run_name
        / "hcc_st/step5_corneto_per_receptor.csv"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", choices=REPEAT_RUN_NAMES)
    arguments = parser.parse_args()
    config = load_config()
    environment = configure_runtime(config)
    if arguments.worker is not None:
        output = run_once(config, arguments.worker)
        print(f"CORNETO worker complete: {arguments.worker} -> {output}")
        return

    validate_signal_extractor(config)
    first_path = run_isolated_worker(config, "run1", environment)
    second_path = run_isolated_worker(config, "run2", environment)
    first = pd.read_csv(first_path)
    second = pd.read_csv(second_path)
    first_diagnostic = pd.read_csv(
        first_path.parent / "step5_corneto_signal_vs_flow_diagnostic.csv"
    )
    second_diagnostic = pd.read_csv(
        second_path.parent / "step5_corneto_signal_vs_flow_diagnostic.csv"
    )
    first_lexicographic = pd.read_csv(
        first_path.parent / "step5_corneto_lexicographic_audit.csv"
    )
    second_lexicographic = pd.read_csv(
        second_path.parent / "step5_corneto_lexicographic_audit.csv"
    )
    first_keys = route_keys(first)
    second_keys = route_keys(second)
    route_exact = first.equals(second)
    lexicographic_exact = first_lexicographic.equals(second_lexicographic)
    solver_status_exact = (
        lexicographic_exact
        and strict_boolean_series(
            first_lexicographic["solver_all_optimal"],
            "run1 solver_all_optimal",
        ).all()
        and strict_boolean_series(
            second_lexicographic["solver_all_optimal"],
            "run2 solver_all_optimal",
        ).all()
        and first_lexicographic["solver_statuses"].astype(str).eq(
            "optimal"
        ).all()
        and second_lexicographic["solver_statuses"].astype(str).eq(
            "optimal"
        ).all()
        and pd.to_numeric(
            first_lexicographic["master_seed"],
            errors="coerce",
        ).eq(config.master_seed).all()
        and pd.to_numeric(
            second_lexicographic["master_seed"],
            errors="coerce",
        ).eq(config.master_seed).all()
        and pd.to_numeric(
            first_lexicographic["solver_seed"],
            errors="coerce",
        ).eq(config.corneto_solver_seed).all()
        and pd.to_numeric(
            second_lexicographic["solver_seed"],
            errors="coerce",
        ).eq(config.corneto_solver_seed).all()
        and pd.to_numeric(
            first_lexicographic["threads"],
            errors="coerce",
        ).eq(config.threads).all()
        and pd.to_numeric(
            second_lexicographic["threads"],
            errors="coerce",
        ).eq(config.threads).all()
    )
    maximum_mip_gap = max(
        float(first_lexicographic["maximum_mip_gap"].max()),
        float(second_lexicographic["maximum_mip_gap"].max()),
    )
    exact = (
        route_exact
        and lexicographic_exact
        and solver_status_exact
        and 0.0 <= maximum_mip_gap <= 1e-9
    )
    repeatability = pd.DataFrame(
        [
            {
                "step": "CORNETO",
                "status": "exact" if exact else "blocked_non_deterministic",
                "run1_edges": len(first),
                "run2_edges": len(second),
                "shared_edge_keys": len(first_keys & second_keys),
                "run1_only_edge_keys": len(first_keys - second_keys),
                "run2_only_edge_keys": len(second_keys - first_keys),
                "run1_flow_only_edges": int(
                    first_diagnostic["flow_only_edges"].sum()
                ),
                "run2_flow_only_edges": int(
                    second_diagnostic["flow_only_edges"].sum()
                ),
                "run1_signal_without_flow_edges": int(
                    first_diagnostic["signal_without_flow_edges"].sum()
                ),
                "run2_signal_without_flow_edges": int(
                    second_diagnostic["signal_without_flow_edges"].sum()
                ),
                "flow_diagnostic_exact": first_diagnostic.equals(
                    second_diagnostic
                ),
                "route_table_exact": route_exact,
                "lexicographic_audit_exact": lexicographic_exact,
                "extraction_basis": (
                    "reaction_sends_activation_c0+"
                    "reaction_sends_inhibition_c0"
                ),
                "flow_used_for_topology": False,
                "secondary_tie_break": (
                    "primary_lock_then_min_signal_edges_then_canonical_first_"
                    "edges_then_activation_first"
                ),
                "solver_status_exact": solver_status_exact,
                "maximum_mip_gap": maximum_mip_gap,
                "master_seed": config.master_seed,
                "solver_seed": config.corneto_solver_seed,
                "threads": config.threads,
            }
        ]
    )
    repeatability.to_csv(
        config.audit_root / "spatial_corneto_repeatability.tsv",
        sep="\t",
        index=False,
    )
    solver_audit = pd.concat(
        [
            first_lexicographic.assign(run="run1"),
            second_lexicographic.assign(run="run2"),
        ],
        ignore_index=True,
    )
    solver_audit["status"] = np.where(
        strict_boolean_series(
            solver_audit["solver_all_optimal"],
            "solver_all_optimal",
        )
        & solver_audit["solver_statuses"].astype(str).eq("optimal")
        & pd.to_numeric(
            solver_audit["maximum_mip_gap"],
            errors="coerce",
        ).between(0.0, 1e-9, inclusive="both")
        & pd.to_numeric(
            solver_audit["master_seed"],
            errors="coerce",
        ).eq(config.master_seed)
        & pd.to_numeric(
            solver_audit["solver_seed"],
            errors="coerce",
        ).eq(config.corneto_solver_seed)
        & pd.to_numeric(
            solver_audit["threads"],
            errors="coerce",
        ).eq(config.threads),
        "PASS",
        "FAIL",
    )
    solver_audit.to_csv(
        config.audit_root / "spatial_corneto_solver_audit.tsv",
        sep="\t",
        index=False,
    )
    (config.corneto_root / "solver_versions.json").write_text(
        json.dumps(version_info(config), indent=2),
        encoding="utf-8",
    )
    for path in (first_path.parent).glob("step5_corneto_*"):
        if path.is_file():
            shutil.copy2(path, config.corneto_root / path.name)
    if not exact or not solver_audit["status"].eq("PASS").all():
        raise RuntimeError(
            "CORNETO lexicographic signal topology changed across independent "
            "fixed-seed repeats."
        )
    print(
        f"CORNETO repeatable lexicographic signal topology: {len(first)} route edges."
    )


if __name__ == "__main__":
    main()
