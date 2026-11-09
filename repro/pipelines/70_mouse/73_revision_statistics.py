#!/usr/bin/env python3
"""Revision analysis of the 35 experimental animals, on the original gram scale.

The original family contains all 21 tumour-weight pairwise comparisons. Holm
adjustment retains that family; the requested triple-versus-doublet contrasts
are extracts, not retrospectively prespecified tests. Terminal body weight is
descriptive only. No baseline weight or longitudinal toxicity endpoint exists.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import logging
import platform
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

import numpy as np
import openpyxl
import scipy
from numpy.typing import NDArray
from scipy import stats

from mouse_population import Population, default_config, load_population


@dataclass(frozen=True)
class AnalysisConfig:
    project_root: Path
    output_dir: Path
    alpha: float = 0.05
    groups: Tuple[str, ...] = tuple("ABCDEFG")
    treatments: Tuple[str, ...] = (
        "Vehicle", "Hydronidone", "Lenvatinib", "Trametinib",
        "Lenvatinib + Trametinib", "Lenvatinib + Hydronidone",
        "Lenvatinib + Trametinib + Hydronidone",
    )
    drug_factors: Tuple[Tuple[int, int, int], ...] = (
        (0, 0, 0), (0, 0, 1), (1, 0, 0), (0, 1, 0),
        (1, 1, 0), (1, 0, 1), (1, 1, 1),
    )


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, records: Sequence[Mapping[str, object]]) -> None:
    if not records:
        raise ValueError(f"Refusing to write an empty analysis table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def holm_adjust(p_values: Sequence[float]) -> NDArray[np.float64]:
    """Strong family-wise error control for the unchanged pairwise test family."""
    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1 or not np.all(np.isfinite(values)):
        raise ValueError("Holm adjustment requires finite one-dimensional p values.")
    if np.any(values < 0) or np.any(values > 1):
        raise ValueError("P values must be between zero and one.")
    order = np.argsort(values, kind="stable")
    corrected = np.minimum(
        1.0, np.maximum.accumulate(values[order] * np.arange(len(values), 0, -1))
    )
    restored = np.empty_like(corrected)
    restored[order] = corrected
    return restored


def linear_contrast(
    values: Mapping[str, NDArray[np.float64]],
    coefficients: Mapping[str, float],
    alpha: float,
    interval_family_size: int = 1,
) -> Dict[str, float]:
    """Independent-group contrast with Welch-Satterthwaite uncertainty.

    Var(contrast) = sum(c_g^2 * s_g^2 / n_g). The approximate degrees of freedom
    are Var(contrast)^2 / sum((c_g^2 * s_g^2 / n_g)^2 / (n_g - 1)).
    Normal sampling is assumed within each arm; n=5 limits precision and the
    reliability of distributional diagnostics. The outcome scale is unchanged.
    """
    if not coefficients or not 0 < alpha < 1 or interval_family_size < 1:
        raise ValueError("Invalid contrast or confidence-interval configuration.")
    estimate = 0.0
    variance = 0.0
    df_denominator = 0.0
    for group, coefficient in coefficients.items():
        group_values = np.asarray(values[group], dtype=float)
        if len(group_values) < 2 or not np.all(np.isfinite(group_values)):
            raise ValueError(f"Contrast group {group} needs finite replicated values.")
        estimate += coefficient * float(np.mean(group_values))
        term = coefficient ** 2 * float(np.var(group_values, ddof=1)) / len(group_values)
        variance += term
        df_denominator += term ** 2 / (len(group_values) - 1)
    if variance <= 0 or df_denominator <= 0:
        raise ValueError("The contrast has no estimable sampling variance.")
    standard_error = float(np.sqrt(variance))
    degrees_freedom = variance ** 2 / df_denominator
    statistic = estimate / standard_error
    critical = float(stats.t.ppf(1 - alpha / 2, degrees_freedom))
    simultaneous = float(stats.t.ppf(
        1 - alpha / (2 * interval_family_size), degrees_freedom
    ))
    return {
        "estimate_g": estimate,
        "standard_error_g": standard_error,
        "satterthwaite_df": degrees_freedom,
        "t_statistic": statistic,
        "p_two_sided_nominal": float(2 * stats.t.sf(abs(statistic), degrees_freedom)),
        "pointwise_95ci_lower_g": estimate - critical * standard_error,
        "pointwise_95ci_upper_g": estimate + critical * standard_error,
        "bonferroni_family95ci_lower_g": estimate - simultaneous * standard_error,
        "bonferroni_family95ci_upper_g": estimate + simultaneous * standard_error,
    }


def endpoint_values(
    population: Population, endpoint: str, groups: Tuple[str, ...]
) -> Dict[str, NDArray[np.float64]]:
    return {
        group: np.array([
            getattr(animal, endpoint) for animal in population.animals
            if animal.group == group
        ], dtype=float)
        for group in groups
    }


def summarize_endpoints(
    population: Population, config: AnalysisConfig
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for endpoint in ("tumour_weight_g", "terminal_body_weight_g"):
        values = endpoint_values(population, endpoint, config.groups)
        for group, treatment in zip(config.groups, config.treatments):
            observations = values[group]
            n = len(observations)
            mean = float(np.mean(observations))
            sd = float(np.std(observations, ddof=1))
            margin = float(stats.t.ppf(1 - config.alpha / 2, n - 1) * sd / np.sqrt(n))
            rows.append({
                "endpoint": endpoint, "group": group, "treatment": treatment,
                "n": n, "mean_g": mean, "sd_g": sd,
                "median_g": float(np.median(observations)),
                "minimum_g": float(np.min(observations)),
                "maximum_g": float(np.max(observations)),
                "mean_pointwise95ci_lower_g": mean - margin,
                "mean_pointwise95ci_upper_g": mean + margin,
                "interpretation": (
                    "Terminal tumour burden; no acquired-resistance endpoint"
                    if endpoint == "tumour_weight_g" else
                    "Terminal measurement only; not weight change or evidence of safety"
                ),
            })
    return rows


def pairwise_results(
    values: Mapping[str, NDArray[np.float64]], config: AnalysisConfig
) -> List[Dict[str, object]]:
    pairs = list(itertools.combinations(config.groups, 2))
    rows: List[Dict[str, object]] = []
    treatments = dict(zip(config.groups, config.treatments))
    for left, right in pairs:
        result = linear_contrast(values, {left: 1, right: -1}, config.alpha, len(pairs))
        rows.append({
            "group_left": left, "group_right": right,
            "treatment_left": treatments[left], "treatment_right": treatments[right],
            "contrast": f"{left} - {right}",
            "n_left": len(values[left]), "n_right": len(values[right]),
            **result,
            "p_holm_all21": 0.0,
            "multiplicity_family": "Original 21 pairwise tumour-weight comparisons",
            "ci_note": "Pointwise CI is unadjusted; Bonferroni CI covers all 21, not Holm inversion",
        })
    adjusted = holm_adjust([float(row["p_two_sided_nominal"]) for row in rows])
    return [{**row, "p_holm_all21": float(p)} for row, p in zip(rows, adjusted)]


def priority_results(rows: Sequence[Mapping[str, object]]) -> List[Dict[str, object]]:
    priority: List[Dict[str, object]] = []
    for row in rows:
        if row["group_left"] not in ("E", "F") or row["group_right"] != "G":
            continue
        reversed_row = dict(row)
        reversed_row.update({
            "group_left": "G", "group_right": row["group_left"],
            "treatment_left": row["treatment_right"],
            "treatment_right": row["treatment_left"],
            "n_left": row["n_right"], "n_right": row["n_left"],
            "contrast": f"G - {row['group_left']}",
            "estimate_g": -float(row["estimate_g"]),
            "t_statistic": -float(row["t_statistic"]),
            "selection_status": "Requested during revision, not originally prespecified",
        })
        for interval in ("pointwise_95ci", "bonferroni_family95ci"):
            reversed_row[f"{interval}_lower_g"] = -float(row[f"{interval}_upper_g"])
            reversed_row[f"{interval}_upper_g"] = -float(row[f"{interval}_lower_g"])
        priority.append(reversed_row)
    if len(priority) != 2:
        raise ValueError("Both requested triple-versus-doublet comparisons are required.")
    return priority


def design_audit(config: AnalysisConfig) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for group, (l_value, t_value, h_value) in zip(config.groups, config.drug_factors):
        rows.append({
            "group": group, "intercept": 1, "L": l_value, "T": t_value, "H": h_value,
            "LT": l_value * t_value, "LH": l_value * h_value,
            "TH": t_value * h_value, "LTH": l_value * t_value * h_value,
        })
    full = np.array([[row[name] for name in list(row)[1:]] for row in rows])
    conditional = np.array([
        [row[name] for name in ("intercept", "T", "H", "TH")]
        for row in rows if row["L"] == 1
    ])
    audit = {
        "full_three_factor_columns": 8,
        "full_three_factor_rank": int(np.linalg.matrix_rank(full)),
        "alias": "TH and LTH columns are identical; the T+H without L arm is absent",
        "conditional_L1_columns": 4,
        "conditional_L1_rank": int(np.linalg.matrix_rank(conditional)),
        "conditional_contrast": "G - E - F + C",
        "interaction_scope": "T by H interaction conditional on L=1, original gram scale",
        "interpretation": "An exploratory, scale-specific interaction; not full three-drug mechanistic synergy",
    }
    if audit["full_three_factor_rank"] != 7 or audit["conditional_L1_rank"] != 4:
        raise ValueError("Unexpected factorial-design identifiability.")
    return rows, audit


def run_analysis(config: AnalysisConfig) -> None:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    population_config = default_config(config.project_root)
    input_paths = (
        population_config.workbook, population_config.exclusion_csv,
        population_config.mirrored_exclusion_csv,
    )
    hashes_before = {str(path): file_hash(path) for path in input_paths}
    population = load_population(population_config)
    treatments = dict(zip(config.groups, config.treatments))
    animal_rows = [
        {**asdict(animal), "treatment": treatments[animal.group]}
        for animal in population.animals
    ]
    write_csv(config.output_dir / "per_animal_35.csv", animal_rows)
    write_csv(config.output_dir / "terminal_endpoint_summary.csv", summarize_endpoints(population, config))
    tumour = endpoint_values(population, "tumour_weight_g", config.groups)
    pairwise = pairwise_results(tumour, config)
    write_csv(config.output_dir / "tumour_pairwise_welch_holm_all21.csv", pairwise)
    write_csv(config.output_dir / "triple_vs_doublets.csv", priority_results(pairwise))
    interaction = linear_contrast(tumour, {"G": 1, "E": -1, "F": -1, "C": 1}, config.alpha)
    interaction_rows = [{
        "contrast": "G - E - F + C", **interaction,
        "expected_G_under_raw_scale_additivity_g": float(
            np.mean(tumour["E"]) + np.mean(tumour["F"]) - np.mean(tumour["C"])
        ),
        "status": "Revision-requested exploratory conditional interaction; nominal p and CI",
        "ci_note": "One contrast; Bonferroni-family interval equals the pointwise interval",
        "direction": "Positive = higher tumour weight than raw-scale additivity predicts",
        "scope": "Conditional on lenvatinib; missing T+H arm prevents full three-factor interaction",
    }]
    write_csv(config.output_dir / "conditional_L1_T_by_H_interaction.csv", interaction_rows)
    design_rows, design = design_audit(config)
    write_csv(config.output_dir / "factorial_design_matrix.csv", design_rows)
    hashes_after = {str(path): file_hash(path) for path in input_paths}
    if hashes_after != hashes_before:
        raise RuntimeError("Input files changed during the animal analysis.")
    receipt = {
        "revision_date": "2026-09-10", "worksheet_records": population.worksheet_record_count,
        "analysis_animals": len(population.animals),
        "external_records_removed_before_endpoint_parsing": [
            {"group": group, "sample_id": sample_id, "worksheet_row": row}
            for group, sample_id, row in population.external_records
        ],
        "exclusion_ledger": list(population.exclusion_records),
        "n_per_arm": {group: len(tumour[group]) for group in config.groups},
        "input_sha256": hashes_before, "inputs_unchanged_after_run": True,
        "script_sha256": {
            name: file_hash(Path(__file__).with_name(name))
            for name in ("73_revision_statistics.py", "mouse_population.py")
        },
        "runtime": {"python": sys.version, "platform": platform.platform(),
                    "numpy": np.__version__, "scipy": scipy.__version__,
                    "openpyxl": openpyxl.__version__},
        "randomness": "No stochastic operation; no random seed influences these statistics",
        "endpoint_scale": "Original tumour weight in grams; no data-driven transform selection",
        "test_family": "All 21 original tumour-weight pairwise tests; two-sided Welch with Holm",
        "confidence_intervals": "Pointwise 95% and conservative Bonferroni 95% family coverage for all 21",
        "small_sample_assumptions": "Independent animals and approximately normal within-arm outcomes; n=5 limits precision",
        "terminal_body_weight": "Descriptive only; no baseline, change, longitudinal tolerability or safety inference",
        "design": design,
    }
    (config.output_dir / "analysis_audit.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    logging.getLogger(__name__).info(
        "Verified %d animals, seven arms; wrote %s", len(population.animals), config.output_dir
    )


def main() -> None:
    root = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path,
        default=root / "repo/results/70_mouse/revision_20260910",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run_analysis(AnalysisConfig(project_root=root, output_dir=args.output_dir))


if __name__ == "__main__":
    main()
