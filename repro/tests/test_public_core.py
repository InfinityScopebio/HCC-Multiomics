from __future__ import annotations

import csv
import hashlib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PRODUCTION_ROOTS = (
    ROOT / "repro/pipelines",
    ROOT / "repro/legacy",
    ROOT / "repro/scripts",
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
    "CoVarNet": "repro/pipelines/20_scrna_cm02/05_rebuild_covarnet_cm02.R",
    "MuSiC": "repro/pipelines/20_scrna_cm02/28_run_music_rnaseq.R",
    "RCTD": "repro/pipelines/50_spatial_cci/10_rctd_from_raw.R",
    "COMMOT": "repro/pipelines/50_spatial_cci/51_commot_from_raw.py",
    "CORNETO": "repro/pipelines/50_spatial_cci/80_corneto_mapk.py",
    "MOFA": "repro/pipelines/60_serum_mofa/15_train_mofa_model.R",
    "pathway_scoring": "repro/pipelines/30_tissue_multiomics/10_progeny.R",
    "custom_projection": "repro/pipelines/20_scrna_cm02/30_build_bulk_cm02_trajectory.R",
}


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PublicCoreTests(unittest.TestCase):
    def production_files(self):
        files = [ROOT / "repro/run_core.sh"]
        for directory in PRODUCTION_ROOTS:
            files.extend(
                path for path in directory.rglob("*")
                if path.is_file() and path.suffix in CODE_SUFFIXES
            )
        return files

    def test_required_method_scripts_exist(self):
        for method, relative in METHOD_SCRIPTS.items():
            with self.subTest(method=method):
                self.assertTrue((ROOT / relative).is_file(), relative)

    def test_production_code_has_no_plotting_calls(self):
        offenders = []
        for path in self.production_files():
            text = path.read_text(encoding="utf-8")
            for token in FORBIDDEN_TOKENS:
                if token in text:
                    offenders.append(f"{path.relative_to(ROOT)}: {token}")
        self.assertEqual([], offenders)

    def test_no_binary_data_or_figure_assets_in_snapshot(self):
        forbidden_suffixes = {
            ".rds", ".h5", ".h5ad", ".hdf5", ".xlsx", ".xls",
            ".pdf", ".png", ".jpg", ".jpeg", ".pptx", ".docx",
        }
        offenders = []
        for path in ROOT.rglob("*"):
            if not path.is_file() or ".git" in path.parts:
                continue
            if "repo/data" in path.as_posix() or "repo/results" in path.as_posix():
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
