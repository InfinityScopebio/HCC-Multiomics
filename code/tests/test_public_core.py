from __future__ import annotations

import ast
import csv
import hashlib
import os
import subprocess
import sys
import shutil
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

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
        self.assertEqual(len(rows), len({row["path"] for row in rows}))
        for row in rows:
            with self.subTest(path=row["path"]):
                path = ROOT / row["path"]
                self.assertTrue(path.is_file())
                self.assertEqual(row["public_sha256"], file_hash(path))
                self.assertEqual(int(row["bytes"]), path.stat().st_size)

    def test_fresh_defaults_and_force_options(self):
        atlas = (ROOT / "code/pipelines/20_scrna_cm02/03_rebuild_atlas_from_local_raw.R").read_text(encoding="utf-8")
        subclusters = (ROOT / "code/pipelines/20_scrna_cm02/04_rebuild_atlas_panel_data.R").read_text(encoding="utf-8")
        runner = (ROOT / "run_core.sh").read_text(encoding="utf-8")
        self.assertIn('resume_enabled("LIUF_RESUME_ATLAS")', atlas)
        self.assertIn('resume_enabled("LIUF_RESUME_SUBCLUSTERS")', subclusters)
        self.assertTrue((ROOT / "code/lib/R/runtime_flags.R").is_file())
        self.assertIn("10_rctd_from_raw.R --force", runner)
        self.assertIn("51_commot_from_raw.py --force", runner)
        self.assertIn("verify_spatial_runtime", runner)
        self.assertLess(
            runner.index("49_cellchat_database.py"),
            runner.index("50_export_scrna_profiles.R"),
        )

    def test_invalid_spatial_interpreter_fails_before_analysis(self):
        env = dict(os.environ)
        env["LIUF_C2L_PYTHON"] = "/nonexistent/spatial-python"
        result = subprocess.run(
            ["bash", "run_core.sh", "spatial"],
            cwd=ROOT,
            env=env,
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not executable", result.stderr)

    def test_commot_cache_provenance_contract(self):
        script = (ROOT / "code/pipelines/50_spatial_cci/51_commot_from_raw.py").read_text(encoding="utf-8")
        for token in (
            "visium_matrix_sha256",
            "positions_sha256",
            "rctd_weights_sha256",
            "cm_recovery_sha256",
            "cellchat_db_sha256",
            "script_sha256",
            "section_seed",
            "commot_version",
            "scanpy_version",
            "cci_input_provenance",
        ):
            with self.subTest(token=token):
                self.assertIn(token, script)
        self.assertIn("cache provenance changed; recomputing", script)
        main_node = next(
            node for node in ast.parse(script).body
            if isinstance(node, ast.FunctionDef) and node.name == "main"
        )

        def run_cache(cache_metadata, expected_provenance, force=False):
            with tempfile.TemporaryDirectory(dir=ROOT) as temporary:
                temporary_root = Path(temporary)
                cache_path = temporary_root / "commot_S1.h5ad"
                cache_path.write_text("test cache", encoding="utf-8")
                config = SimpleNamespace(
                    tumor_sections=("S1",),
                    commot_cache=lambda section: cache_path,
                    commot_root=temporary_root,
                    result_root=temporary_root,
                    audit_root=temporary_root,
                )
                calls = []

                def fake_run_section(config, section, output_root, database):
                    calls.append(section)
                    return cache_path, {"section": section, "n_spots": 30}

                def write_table(rows, path, sep, index):
                    Path(path).write_text(repr(rows), encoding="utf-8")

                namespace = {
                    "argparse": __import__("argparse"),
                    "load_config": lambda: config,
                    "configure_runtime": lambda config: None,
                    "cellchat_database": lambda config: object(),
                    "input_provenance": lambda config, section: expected_provenance,
                    "sc": SimpleNamespace(
                        read_h5ad=lambda path: SimpleNamespace(uns=cache_metadata)
                    ),
                    "run_section": fake_run_section,
                    "numeric_fingerprint": lambda config, cached, section, cache_path: {
                        "section": section,
                        "n_spots": 30,
                    },
                    "pd": SimpleNamespace(
                        DataFrame=lambda rows: SimpleNamespace(
                            to_csv=lambda path, sep, index: write_table(rows, path, sep, index)
                        )
                    ),
                }
                exec(compile(ast.Module(body=[main_node], type_ignores=[]), "commot_main", "exec"), namespace)
                argv = ["51_commot_from_raw.py"] + (["--force"] if force else [])
                with mock.patch.object(sys, "argv", argv):
                    namespace["main"]()
                return calls

        matching = {"cci_input_provenance": {"input": "v1"}}
        self.assertEqual([], run_cache(matching, {"input": "v1"}))
        self.assertEqual(["S1"], run_cache({}, {"input": "v1"}))
        self.assertEqual(["S1"], run_cache(matching, {"input": "v2"}))
        self.assertEqual(["S1"], run_cache(matching, {"input": "v1"}, force=True))

    def test_raw_mofa_dispatch_and_provenance(self):
        runner = (ROOT / "code/pipelines/60_serum_mofa/run_all.sh").read_text(encoding="utf-8")
        self.assertIn('if [[ "${MOFA_ANALYSIS_SOURCE}" != "raw_derived" ]]; then\n  run_stage "05_materialize_historical_boundary.R"\nfi', runner)
        self.assertIn('if [[ "${MOFA_ANALYSIS_SOURCE}" != "raw_derived" ]]; then\n  run_stage "18_article_k_scan_and_variance.R"\nfi', runner)
        self.assertIn('if [[ "${MOFA_ANALYSIS_SOURCE}" == "raw_derived" ]]; then\n  run_stage "61_finalize_raw_audit.R"', runner)
        raw = (ROOT / "code/pipelines/60_serum_mofa/61_finalize_raw_audit.R").read_text(encoding="utf-8")
        self.assertIn("require_analysis_mofa_boundary(cfg)", raw)
        self.assertIn("ARTICLE_EQUIVALENCE_NOT_ESTABLISHED", raw)
        self.assertIn("mofa_raw_verdict.tsv", raw)
        self.assertIn('fast1000_replay = "k_scan_fast1000_replay"', raw)
        self.assertIn('slow2000_replay = "final_slow2000_replay"', raw)
        self.assertIn("profile = unname(metadata_profiles[[profile]])", raw)
        mofa_config = (ROOT / "code/pipelines/60_serum_mofa/config.R").read_text(encoding="utf-8")
        self.assertIn('model_dir <- file.path(result_dir, "models")', mofa_config)
        worker = (ROOT / "code/pipelines/60_serum_mofa/16_train_fast_worker.sh").read_text(encoding="utf-8")
        self.assertIn('RESULT_DIR="${CANONICAL_RESULT_DIR}/raw_to_figure"', worker)
        for script_name in (
            "25_rebuild_tcga_axes.R",
            "30_tcga_survival.R",
            "50_supplemental_tables.R",
        ):
            text = (ROOT / f"code/pipelines/60_serum_mofa/{script_name}").read_text(encoding="utf-8")
            with self.subTest(script=script_name):
                self.assertIn("source_model_sha256", text)
                self.assertIn("source_input_sha256", text)

        with tempfile.TemporaryDirectory(dir=ROOT) as temporary:
            sandbox = Path(temporary)
            shutil.copytree(ROOT / "code", sandbox / "code")
            shutil.copytree(ROOT / "config", sandbox / "config")
            spy = sandbox / "rscript-spy"
            spy.write_text(
                "#!/usr/bin/env bash\n"
                "set -euo pipefail\n"
                "if [[ \"$*\" == *LIUF_R_RUNTIME_OK* ]]; then\n"
                "  printf 'LIUF_R_RUNTIME_OK'\n"
                "  exit 0\n"
                "fi\n"
                "basename \"$1\" >>\"${SPY_LOG}\"\n",
                encoding="utf-8",
            )
            spy.chmod(0o755)
            for mode in ("article", "raw_derived"):
                with self.subTest(mode=mode):
                    log = sandbox / f"{mode}.calls"
                    env = dict(os.environ)
                    env.update({
                        "LIUF_RSCRIPT": str(spy),
                        "LIUF_PYTHON": sys.executable,
                        "SPY_LOG": str(log),
                        "MOFA_ANALYSIS_SOURCE": mode,
                        "MOFA_RUN_RAW_DIAGNOSTIC": "0",
                        "MOFA_FORCE_RETRAIN": "0",
                        "MOFA_WORKERS": "1",
                    })
                    result = subprocess.run(
                        ["bash", "code/pipelines/60_serum_mofa/run_all.sh"],
                        cwd=sandbox,
                        env=env,
                        capture_output=True,
                        text=True,
                        timeout=30,
                        check=False,
                    )
                    self.assertEqual(result.returncode, 0, result.stderr)
                    calls = log.read_text(encoding="utf-8").splitlines()
                    if mode == "raw_derived":
                        self.assertEqual(calls, [
                            "10_rebuild_incdec_input.R",
                            "25_rebuild_tcga_axes.R",
                            "30_tcga_survival.R",
                            "40_factor2_classifier.R",
                            "50_supplemental_tables.R",
                            "61_finalize_raw_audit.R",
                        ])
                        runtime = sandbox / "results/60_serum_mofa/raw_to_figure/pipeline_stage_runtime.tsv"
                    else:
                        self.assertEqual(calls, [
                            "05_materialize_historical_boundary.R",
                            "10_rebuild_incdec_input.R",
                            "18_article_k_scan_and_variance.R",
                            "25_rebuild_tcga_axes.R",
                            "30_tcga_survival.R",
                            "40_factor2_classifier.R",
                            "50_supplemental_tables.R",
                            "61_finalize_core_audit.R",
                        ])
                        runtime = sandbox / "results/60_serum_mofa/pipeline_stage_runtime.tsv"
                    self.assertTrue(runtime.is_file())

            env = dict(os.environ)
            env.update({
                "LIUF_RSCRIPT": str(spy),
                "LIUF_PYTHON": sys.executable,
                "SPY_LOG": str(sandbox / "invalid.calls"),
                "MOFA_ANALYSIS_SOURCE": "article",
                "MOFA_RUN_RAW_DIAGNOSTIC": "1",
                "MOFA_FORCE_RETRAIN": "1",
                "MOFA_WORKERS": "1",
            })
            invalid = subprocess.run(
                ["bash", "code/pipelines/60_serum_mofa/run_all.sh"],
                cwd=sandbox,
                env=env,
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
            self.assertNotEqual(invalid.returncode, 0)
            self.assertIn("requires MOFA_ANALYSIS_SOURCE=raw_derived", invalid.stderr)
            self.assertFalse((sandbox / "invalid.calls").exists())

    def test_raw_audit_contract_executes(self):
        fixture = ROOT / "code/tests/test_raw_audit_contract.R"
        self.assertTrue(fixture.is_file())
        text = fixture.read_text(encoding="utf-8")
        for token in (
            "diagnostic_difference",
            "repeatability_failure",
            "another_model",
            "duplicated_ids",
            "missing_id",
            "file_hash_difference",
            "semantic_difference",
            "nonfinite",
            "empty",
            "wrong_downstream_id",
            "assert_article_audit_contract",
            "empty_article_input",
            "empty_article_boundary",
        ):
            with self.subTest(token=token):
                self.assertIn(token, text)
        result = subprocess.run(
            [os.environ.get("LIUF_RSCRIPT", "Rscript"), str(fixture)],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Raw and article audit contract fixtures PASS", result.stdout)

    def test_storm_release_is_pinned(self):
        script = (ROOT / "code/pipelines/40_tissue_resistance/41_rebuild_storm_gsea.R").read_text(encoding="utf-8")
        self.assertIn('STORM_MSIGDB_VERSION", unset = "2025.1.Hs"', script)
        self.assertIn("observed_msigdb_version, storm_msigdb_version", script)
        self.assertIn("storm_msigdb_cache", script)
        self.assertIn("no download is permitted", script)
        self.assertIn("observed_msigdb_version", script)
        self.assertIn("storm_msigdb_version.tsv", script)
        mofa = (ROOT / "code/pipelines/60_serum_mofa/25_rebuild_tcga_axes.R").read_text(encoding="utf-8")
        self.assertIn('MOFA_MSIGDB_VERSION", unset = "2025.1.Hs"', mofa)
        self.assertIn("observed_version, msigdb_version", mofa)
        self.assertIn("msigdb_cache", mofa)
        self.assertIn("no download is permitted", mofa)
        self.assertIn("mofa_msigdb_version.tsv", mofa)

    def test_manifest_inclusion_contract(self):
        with (ROOT / "config/manifests/panel_registry.tsv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(rows)
        included = []
        for row in rows:
            for relative in row["included_generator_paths"].split(";"):
                if relative:
                    included.append(relative)
                    with self.subTest(path=relative):
                        self.assertTrue((ROOT / relative).is_file(), relative)
            self.assertEqual(row["accepted_artifact_status"], "historical_reference_only")
        self.assertTrue(included)
        for relative in (
            "config/manifests/cm02_historical_downstream.tsv",
            "config/manifests/supplemental2_rank_survey_historical.tsv",
            "code/pipelines/20_scrna_cm02/10_materialize_historical_downstream.R",
            "code/pipelines/20_scrna_cm02/20_materialize_supplemental2_rank.R",
        ):
            with self.subTest(path=relative):
                self.assertTrue((ROOT / relative).is_file(), relative)

    def test_license_and_stale_configuration(self):
        self.assertTrue((ROOT / "LICENSE").is_file())
        self.assertTrue((ROOT / "NOTICE.md").is_file())
        boundaries = (ROOT / "docs/BOUNDARIES.md").read_text(encoding="utf-8")
        self.assertNotIn("re" + "po/audit", boundaries)
        spatial_config = (ROOT / "code/pipelines/50_spatial_cci/config.py").read_text(encoding="utf-8")
        self.assertNotIn("historical_plot_seed", spatial_config)
        self.assertNotIn("panel_root", spatial_config)
        mofa_config = (ROOT / "code/pipelines/60_serum_mofa/config.R").read_text(encoding="utf-8")
        self.assertNotIn("panel_dir", mofa_config)
        project_config = (ROOT / "code/lib/python/project_config.py").read_text(encoding="utf-8")
        self.assertNotIn("render_error_fraction", project_config)

    def test_input_ignore_patterns_cover_symlinks(self):
        ignore = (ROOT / ".gitignore").read_text(encoding="utf-8")
        for pattern in (
            "data/private",
            "data/public",
            "data/validation",
            "data/curation/*",
            "data/reference/*",
            "!data/curation/resistance/",
            "!data/reference/annotations/",
            "!data/reference/KEGG/",
        ):
            with self.subTest(pattern=pattern):
                self.assertIn(pattern, ignore.splitlines())

    def test_environment_documentation(self):
        environment = (ROOT / "code/env/README.md").read_text(encoding="utf-8")
        self.assertIn("Python 3.11 or newer", environment)
        self.assertIn("R 4.5.2", environment)
        self.assertIn("scanpy 1.11.5", environment)
        self.assertIn("COMMOT 0.0.3", environment)
        self.assertIn("CORNETO 0.9.1a6", environment)
        self.assertTrue((ROOT / "code/env/R_libpaths.txt").is_file())

    def test_python_syntax_without_bytecode(self):
        for path in self.production_files():
            if path.suffix != ".py":
                continue
            with self.subTest(path=path.relative_to(ROOT)):
                compile(path.read_text(encoding="utf-8"), str(path), "exec")


if __name__ == "__main__":
    unittest.main()
