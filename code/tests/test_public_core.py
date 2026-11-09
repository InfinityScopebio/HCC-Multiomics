from __future__ import annotations

import csv
import hashlib
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PRODUCTION_ROOTS = (
    ROOT / "code/pipelines",
    ROOT / "code/vendor",
    ROOT / "code/scripts",
    ROOT / "code/lib",
)
CODE_SUFFIXES = {".R", ".r", ".py", ".sh"}
FORBIDDEN_TOKENS = (
    "ggplot2",
    "ComplexHeatmap",
    "enrichplot",
    "survminer",
    "cowplot",
    "patchwork",
    "ggsave",
    "cairo_pdf",
    "matplotlib",
    "seaborn",
    "plotly",
    "savefig",
    "plot_pairwise_network",
    "plot_river",
)
METHOD_SCRIPTS = {
    "CoVarNet": "code/pipelines/20_scrna_cm02/05_rebuild_covarnet_cm02.R",
    "MuSiC": "code/pipelines/20_scrna_cm02/28_run_music_rnaseq.R",
    "RCTD": "code/pipelines/50_spatial_cci/10_rctd_from_raw.R",
    "COMMOT": "code/pipelines/50_spatial_cci/51_commot_from_raw.py",
    "CORNETO": "code/pipelines/50_spatial_cci/80_corneto_mapk.py",
    "MOFA": "code/pipelines/60_serum_mofa/15_train_mofa_model.R",
    "pathway_scoring": "code/pipelines/30_tissue_multiomics/10_progeny.R",
    "custom_projection": "code/pipelines/20_scrna_cm02/30_build_bulk_cm02_trajectory.R",
}
EXPECTED_ROOT_ENTRIES = {
    "README.md",
    "run_core.sh",
    "code",
    "config",
    "docs",
    "data",
    "SOURCE_MANIFEST.tsv",
}


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PublicCoreTests(unittest.TestCase):
    def production_files(self):
        files = [ROOT / "run_core.sh"]
        for directory in PRODUCTION_ROOTS:
            files.extend(
                path for path in directory.rglob("*")
                if path.is_file() and path.suffix in CODE_SUFFIXES
            )
        return files

    def test_unified_root_layout(self):
        for name in EXPECTED_ROOT_ENTRIES:
            with self.subTest(path=name):
                self.assertTrue((ROOT / name).exists(), name)
        for name in ("re" + "pro", "re" + "po", "mofa_v5_mefisto", "plan"):
            with self.subTest(removed_root=name):
                self.assertFalse((ROOT / name).exists(), name)

    def test_readme_has_no_scope_labels(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        for phrase in (
            "Computation-only code snapshot",
            "This repository intentionally contains",
            "Reviewer 2, minor comment",
        ):
            with self.subTest(phrase=phrase):
                self.assertNotIn(phrase, readme)

    def test_required_method_scripts_exist(self):
        for method, relative in METHOD_SCRIPTS.items():
            with self.subTest(method=method):
                self.assertTrue((ROOT / relative).is_file(), relative)
        self.assertFalse((ROOT / "code/pipelines/70_mouse").exists())

    def test_configuration_and_key_paths(self):
        sys.path.insert(0, str(ROOT / "code/lib/python"))
        from project_config import load_project_config

        config = load_project_config()
        self.assertEqual(ROOT, config.paths.root)
        self.assertEqual(ROOT, config.paths.repo)
        self.assertEqual(ROOT / "code", config.paths.repro)
        self.assertEqual(ROOT / "data", config.paths.data)
        self.assertEqual(ROOT / "results", config.paths.results)
        self.assertEqual(ROOT / "audit", config.paths.audit)
        self.assertEqual(ROOT / "logs", config.paths.logs)
        self.assertEqual(ROOT / "manifests", config.paths.manifests)
        mofa_config = (ROOT / "code/pipelines/60_serum_mofa/config.R").read_text(
            encoding="utf-8"
        )
        self.assertIn('"data",\n  "reference",\n  "KEGG",', mofa_config)
        self.assertNotIn('"mofa_v5_mefisto",', mofa_config)
        tissue_preflight = (
            ROOT / "code/pipelines/30_tissue_multiomics/00_validate_inputs.R"
        ).read_text(encoding="utf-8")
        self.assertIn('"config", "manifests", "panel_registry.tsv"', tissue_preflight)

    def test_no_stale_runtime_paths(self):
        offenders = []
        for path in self.production_files():
            if "tests" in path.parts:
                continue
            text = path.read_text(encoding="utf-8")
            for token in ("repro/", "repo/"):
                if token in text:
                    offenders.append(f"{path.relative_to(ROOT)}: {token}")
        self.assertEqual([], offenders)

    def test_production_code_has_no_plotting_calls(self):
        offenders = []
        for path in self.production_files():
            text = path.read_text(encoding="utf-8")
            for token in FORBIDDEN_TOKENS:
                if token in text:
                    offenders.append(f"{path.relative_to(ROOT)}: {token}")
        self.assertEqual([], offenders)

    def test_no_binary_data_or_figure_assets_in_repository(self):
        forbidden_suffixes = {
            ".rds", ".h5", ".h5ad", ".hdf5", ".xlsx", ".xls",
            ".pdf", ".png", ".jpg", ".jpeg", ".pptx", ".docx",
        }
        offenders = []
        for path in ROOT.rglob("*"):
            if not path.is_file() or ".git" in path.parts:
                continue
            if "data" in path.parts or "results" in path.parts:
                continue
            if path.suffix.lower() in forbidden_suffixes:
                offenders.append(str(path.relative_to(ROOT)))
        self.assertEqual([], offenders)

    def test_source_manifest_hashes(self):
        manifest = ROOT / "SOURCE_MANIFEST.tsv"
        self.assertTrue(manifest.is_file())
        with manifest.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertGreater(len(rows), 50)
        for row in rows:
            with self.subTest(path=row["path"]):
                path = ROOT / row["path"]
                self.assertTrue(path.is_file())
                self.assertEqual(row["public_sha256"], file_hash(path))

    def test_python_syntax_without_bytecode(self):
        for path in self.production_files():
            if path.suffix != ".py":
                continue
            with self.subTest(path=path.relative_to(ROOT)):
                compile(path.read_text(encoding="utf-8"), str(path), "exec")


if __name__ == "__main__":
    unittest.main()
