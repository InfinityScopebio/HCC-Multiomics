"""Verify exclusion behavior and statistics against independent implementations."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from typing import Any, List

import numpy as np
import openpyxl
from scipy import stats
from statsmodels.stats.multitest import multipletests

from mouse_population import default_config, load_population


def load_statistics() -> Any:
    path = Path(__file__).with_name("73_revision_statistics.py")
    spec = importlib.util.spec_from_file_location("mouse_revision_statistics", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot import mouse statistics module.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class AnimalRevisionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.module = load_statistics()
        self.root = Path(__file__).resolve().parents[3]
        self.population_config = default_config(self.root)
        self.population = load_population(self.population_config)

    def test_population_excludes_nonexperimental_record_from_both_endpoints(self) -> None:
        self.assertEqual(len(self.population.animals), 35)
        self.assertEqual(self.population.worksheet_record_count, 36)
        self.assertNotIn(21, [animal.sample_id for animal in self.population.animals])
        for endpoint in ("tumour_weight_g", "terminal_body_weight_g"):
            values = self.module.endpoint_values(self.population, endpoint, tuple("ABCDEFG"))
            self.assertEqual([len(values[group]) for group in "ABCDEFG"], [5] * 7)
        self.assertAlmostEqual(np.mean([
            animal.tumour_weight_g for animal in self.population.animals if animal.group == "E"
        ]), 2.82)

    def test_external_endpoint_values_are_not_parsed(self) -> None:
        with tempfile.TemporaryDirectory(prefix="mouse_revision_") as directory:
            path = Path(directory) / "external_values_unreadable.xlsx"
            workbook = openpyxl.load_workbook(self.population_config.workbook)
            sheet = workbook.worksheets[0]
            for row in sheet.iter_rows(min_row=2):
                if row[3].value == 21:
                    row[4].value = "not an endpoint from this experiment"
                    row[5].value = None
            workbook.save(path)
            workbook.close()
            actual = load_population(replace(self.population_config, workbook=path))
            self.assertEqual(actual.animals, self.population.animals)

    def test_changed_or_inconsistent_exclusion_ledger_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="mouse_revision_") as directory:
            path = Path(directory) / "ledger.csv"
            path.write_text(
                self.population_config.exclusion_csv.read_text().replace(
                    "2026-09-10", "2026-09-09"
                )
            )
            with self.assertRaisesRegex(ValueError, "ledgers differ"):
                load_population(replace(self.population_config, exclusion_csv=path))

    def test_pairwise_contrasts_match_scipy_welch_and_holm_matches_statsmodels(self) -> None:
        config = self.module.AnalysisConfig(self.root, self.root)
        values = self.module.endpoint_values(self.population, "tumour_weight_g", config.groups)
        results = self.module.pairwise_results(values, config)
        self.assertEqual(len(results), 21)
        independent_p: List[float] = []
        for result in results:
            test = stats.ttest_ind(
                values[result["group_left"]], values[result["group_right"]], equal_var=False
            )
            independent_p.append(float(test.pvalue))
            self.assertAlmostEqual(result["p_two_sided_nominal"], test.pvalue, places=13)
            self.assertAlmostEqual(result["satterthwaite_df"], test.df, places=12)
            interval = test.confidence_interval(confidence_level=0.95)
            self.assertAlmostEqual(result["pointwise_95ci_lower_g"], interval.low, places=12)
            self.assertAlmostEqual(result["pointwise_95ci_upper_g"], interval.high, places=12)
        expected = multipletests(independent_p, method="holm")[1]
        np.testing.assert_allclose([row["p_holm_all21"] for row in results], expected, rtol=1e-12)

    def test_interaction_direction_and_identifiability(self) -> None:
        config = self.module.AnalysisConfig(self.root, self.root)
        values = self.module.endpoint_values(self.population, "tumour_weight_g", config.groups)
        result = self.module.linear_contrast(
            values, {"G": 1, "E": -1, "F": -1, "C": 1}, config.alpha
        )
        self.assertAlmostEqual(result["estimate_g"], 2.98)
        self.assertGreater(result["pointwise_95ci_lower_g"], 0)
        _, design = self.module.design_audit(config)
        self.assertEqual(design["full_three_factor_rank"], 7)
        self.assertEqual(design["conditional_L1_rank"], 4)

    def test_priority_comparisons_retain_original_family(self) -> None:
        config = self.module.AnalysisConfig(self.root, self.root)
        values = self.module.endpoint_values(self.population, "tumour_weight_g", config.groups)
        pairwise = self.module.pairwise_results(values, config)
        priority = self.module.priority_results(pairwise)
        self.assertEqual([row["contrast"] for row in priority], ["G - E", "G - F"])
        for row in priority:
            self.assertGreater(row["p_holm_all21"], 0.05)
            self.assertLess(row["bonferroni_family95ci_lower_g"], 0)
            self.assertGreater(row["bonferroni_family95ci_upper_g"], 0)
            self.assertIn("not originally prespecified", row["selection_status"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
