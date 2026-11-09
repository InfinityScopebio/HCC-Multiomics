#!/usr/bin/env python3
"""Immutable paths and scientific scope for the spatial/CCI reproduction."""

from __future__ import annotations

import hashlib
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path


PYTHON_LIB_ROOT = Path(__file__).resolve().parents[2] / "lib" / "python"
sys.path.insert(0, str(PYTHON_LIB_ROOT))

try:
    from project_config import R_SEED_MAX, derive_r_seed
except ModuleNotFoundError as error:
    if error.name != "tomllib":
        raise
    # The c2l Python 3.10 runtime lacks tomllib; this module only needs the
    # shared R seed bound and does not parse the project TOML configuration.
    R_SEED_MAX = 2**31 - 1
    R_SEED_MODULUS = 2**31

    def derive_r_seed(base_seed: int, offset: int = 0) -> int:
        """Derive an R-compatible seed using the shared project contract."""

        if isinstance(base_seed, bool) or not isinstance(base_seed, int):
            raise TypeError("base_seed must be an integer")
        if isinstance(offset, bool) or not isinstance(offset, int):
            raise TypeError("offset must be an integer")
        if not 0 <= base_seed <= R_SEED_MAX:
            raise ValueError(f"base_seed must be in [0, {R_SEED_MAX}]")
        return (base_seed + offset) % R_SEED_MODULUS


@dataclass(frozen=True)
class SpatialConfig:
    """Configuration for article-equivalence and diagnostic spatial runs."""

    project_root: Path
    run_variant: str = "article_equivalent"
    master_seed: int = 20260807
    historical_plot_seed: int = 1
    corneto_solver_seed: int = 0
    threads: int = 1
    commot_cache_version: str = "historical_rctd_cm02_barcode_v4"

    @property
    def rctd_source(self) -> str:
        if self.run_variant == "article_equivalent":
            return "historical_frozen"
        if self.run_variant == "raw_rctd_diagnostic":
            return "raw_recomputed"
        raise ValueError(f"Unsupported spatial run variant: {self.run_variant}")

    @property
    def article_equivalent(self) -> bool:
        return self.run_variant == "article_equivalent"

    @property
    def active_commot_cache_version(self) -> str:
        if self.article_equivalent:
            return self.commot_cache_version
        return "raw_rctd_diagnostic_v1"

    @property
    def pipeline_root(self) -> Path:
        return self.project_root / "code/pipelines/50_spatial_cci"

    @property
    def legacy_code_root(self) -> Path:
        return self.project_root / "code/vendor/source/hcc_st"

    @property
    def raw_spatial_root(self) -> Path:
        return self.project_root / "data/public/raw/HCC-ST"

    @property
    def scrna_result_root(self) -> Path:
        return self.project_root / "results/20_scrna_cm02"

    @property
    def scrna_curation_root(self) -> Path:
        return self.project_root / "data/curation/scRNA"

    @property
    def scrna_atlas(self) -> Path:
        return self.scrna_curation_root / "combined_annotated_v5.rds"

    @property
    def scrna_hepatocyte(self) -> Path:
        return self.scrna_curation_root / "hep_annotated.rds"

    @property
    def historical_rctd_reference(self) -> Path:
        return (
            self.project_root
            / "data/curation/spatial/RCTD/ST3_rctd_reference.rds"
        )

    @property
    def historical_rctd_reference_sha256(self) -> str:
        return (
            "1dd7de4e9d286c89dc22ccc801f31a1bb18e3c090b577d0c25585b93e337eab9"
        )

    @property
    def historical_rctd_sigma_pilot_manifest(self) -> Path:
        return (
            self.project_root
            / "config/manifests/spatial_rctd_historical_sigma_pilot.tsv"
        )

    @property
    def historical_intermediate_manifest(self) -> Path:
        return (
            self.project_root
            / "config/manifests/spatial_historical_intermediates.tsv"
        )

    @property
    def article_cci_validation_manifest(self) -> Path:
        return (
            self.project_root
            / "config/manifests/spatial_article_cci_validation.tsv"
        )

    @property
    def historical_corneto_manifest(self) -> Path:
        return (
            self.project_root
            / "config/manifests/spatial_historical_corneto.tsv"
        )

    @property
    def historical_corneto_root(self) -> Path:
        return self.project_root / "data/curation/spatial/CORNETO"

    @property
    def historical_corneto_route(self) -> Path:
        return self.historical_corneto_root / "step5_corneto_per_receptor.csv"

    @property
    def historical_corneto_measurements(self) -> Path:
        return self.historical_corneto_root / "step5_corneto_measurements.csv"

    @property
    def historical_corneto_validation(self) -> Path:
        return self.audit_root / "spatial_historical_corneto_validation.tsv"

    @property
    def historical_intermediate_validation(self) -> Path:
        return (
            self.audit_root
            / "spatial_historical_intermediate_validation.tsv"
        )

    @property
    def historical_rctd_weight_root(self) -> Path:
        return self.project_root / "data/curation/spatial/RCTD/weights"

    @property
    def historical_covarnet_model(self) -> Path:
        return (
            self.project_root
            / "data/curation/spatial/CoVarNet/cm_network.rds"
        )

    @property
    def msigdb_archive(self) -> Path:
        return (
            self.project_root
            / "data/reference/MSigDB"
            / "msigdb_v2023.2.Hs_files_to_download_locally.zip"
        )

    @property
    def msigdb_archive_sha256(self) -> str:
        return (
            "0aa92b11ad56204905ced89581a65affe7ce422c11c9660851c2ec6829a1c33a"
        )

    @property
    def scrna_curation_audit(self) -> Path:
        return self.project_root / "audit/scrna_frozen_curation_objects.tsv"

    @property
    def scrna_count_audit(self) -> Path:
        return self.project_root / "audit/scrna_full_count_rebuild.tsv"

    @property
    def covarnet_model(self) -> Path:
        return self.scrna_result_root / "covarnet/cm_network_frozen_curated.rds"

    @property
    def covarnet_vendor(self) -> Path:
        return (
            self.project_root
            / "code/pipelines/20_scrna_cm02/covarnet_v030_vendor.R"
        )

    @property
    def spatial_covarnet_root(self) -> Path:
        return self.result_root / "covarnet_reference"

    @property
    def spatial_covarnet_model(self) -> Path:
        return (
            self.spatial_covarnet_root
            / "cm_network_frozen_curated_reference_historical_filter.rds"
        )

    @property
    def covarnet_nmf_audit(self) -> Path:
        return self.project_root / "audit/scrna_covarnet_nmf_rebuild.tsv"

    @property
    def cm02_scope_audit(self) -> Path:
        return self.project_root / "audit/scrna_cm02_scope.tsv"

    @property
    def validation_root(self) -> Path:
        return (
            self.project_root
            / "data/validation/frozen_targets/spatial_cci/hcc_st"
        )

    @property
    def scrna_validation_root(self) -> Path:
        return (
            self.project_root
            / "data/validation/frozen_targets/scRNA"
        )

    @property
    def forbidden_reference_root(self) -> Path:
        return self.project_root / "data/reference/spatial_cci"

    @property
    def result_root(self) -> Path:
        root = self.project_root / "results/50_spatial_cci"
        if self.run_variant == "article_equivalent":
            return root
        if self.run_variant == "raw_rctd_diagnostic":
            return root / "raw_rctd_diagnostic"
        raise ValueError(f"Unsupported spatial run variant: {self.run_variant}")

    @property
    def audit_root(self) -> Path:
        if self.article_equivalent:
            return self.project_root / "audit"
        return self.result_root / "audit"

    @property
    def log_root(self) -> Path:
        return self.result_root / "logs"

    @property
    def rctd_root(self) -> Path:
        return self.result_root / "rctd"

    @property
    def rctd_weight_root(self) -> Path:
        return self.rctd_root / "weights"

    @property
    def rctd_reference(self) -> Path:
        return self.rctd_root / "reference/rctd_reference.rds"

    @property
    def cm_recovery_root(self) -> Path:
        return self.result_root / "cm_recovery"

    @property
    def cm_recovery_file(self) -> Path:
        return self.cm_recovery_root / "CM_recover.rds"

    @property
    def cm_recovery_with_barcode_file(self) -> Path:
        return self.cm_recovery_root / "CM_recover_with_barcode.csv.gz"

    @property
    def cm_recovery_with_barcode_rds(self) -> Path:
        return self.cm_recovery_root / "CM_recover_with_barcode.rds"

    @property
    def radial_root(self) -> Path:
        return self.result_root / "radial"

    @property
    def pathway_root(self) -> Path:
        return self.result_root / "pathways"

    @property
    def pathway_score_file(self) -> Path:
        return self.pathway_root / "pathway_scores_raw_visium.rds"

    @property
    def pathway_merged_file(self) -> Path:
        return self.pathway_root / "pathway_scores_with_cm02.rds"

    @property
    def commot_root(self) -> Path:
        return self.result_root / "commot"

    @property
    def cci_root(self) -> Path:
        return self.result_root / "cci"

    @property
    def corneto_root(self) -> Path:
        return self.result_root / "corneto"

    @property
    def panel_root(self) -> Path:
        return self.result_root / "panels"

    @property
    def validation_result_root(self) -> Path:
        return self.result_root / "validation"

    @property
    def rscript(self) -> Path:
        return Path(
            os.environ.get("LIUF_RSCRIPT", "/usr/local/bin/Rscript")
        )

    @property
    def c2l_python(self) -> Path:
        return Path(
            os.environ.get(
                "LIUF_C2L_PYTHON",
                "/Users/pyramidsnail/miniconda3/envs/c2l/bin/python",
            )
        )

    @property
    def all_sections(self) -> tuple[str, ...]:
        return (
            "HCC-1N",
            "HCC-1L",
            "HCC-1T",
            "HCC-2N",
            "HCC-2L",
            "HCC-2T",
            "HCC-3N",
            "HCC-3L",
            "HCC-3T",
            "HCC-4N",
            "HCC-4L",
            "HCC-4T",
            "HCC-5A",
            "HCC-5B",
            "HCC-5C",
            "HCC-5D",
        )

    @property
    def tumor_sections(self) -> tuple[str, ...]:
        return (
            "HCC-1T",
            "HCC-2T",
            "HCC-3T",
            "HCC-4T",
            "HCC-5A",
            "HCC-5B",
            "HCC-5C",
            "HCC-5D",
        )

    @property
    def cm02_types(self) -> tuple[str, ...]:
        return (
            "HSC",
            "Pericyte",
            "Treg",
            "IL1B_Mp",
            "SPP1_Mp",
            "TREM2_Mp",
            "CXCL10_Mp",
            "Hep_C4",
        )

    @property
    def cm02_covarnet_types(self) -> tuple[str, ...]:
        return (
            "HSC",
            "Pericyte",
            "Treg",
            "IL1B+M_",
            "SPP1+M_",
            "TREM2+M_",
            "CXCL10+M_",
            "Hep_C4",
        )

    def visium_dir(self, section: str) -> Path:
        return self.raw_spatial_root / section

    def weight_rds(self, section: str) -> Path:
        return self.rctd_weight_root / f"{section}_weights.rds"

    def weight_csv(self, section: str) -> Path:
        return self.rctd_weight_root / f"{section}_weights.csv.gz"

    def historical_weight_rds(self, section: str) -> Path:
        if section not in self.all_sections:
            raise ValueError(f"Unknown spatial section: {section}")
        if section.startswith("HCC-5"):
            prefix = "ST3"
        elif section.endswith("T"):
            prefix = "ST5"
        else:
            prefix = "ST9k"
        return self.historical_rctd_weight_root / f"{prefix}_{section}_weights.rds"

    def active_weight_rds(self, section: str) -> Path:
        if self.rctd_source == "raw_recomputed":
            return self.weight_rds(section)
        return self.historical_weight_rds(section)

    def commot_cache(self, section: str) -> Path:
        return self.commot_root / f"commot_{section}.h5ad"

    def ensure_output_dirs(self) -> None:
        for path in (
            self.result_root,
            self.audit_root,
            self.log_root,
            self.rctd_weight_root,
            self.rctd_reference.parent,
            self.cm_recovery_root,
            self.radial_root,
            self.pathway_root,
            self.commot_root,
            self.cci_root,
            self.corneto_root,
            self.panel_root,
            self.validation_result_root,
        ):
            path.mkdir(parents=True, exist_ok=True)

    def assert_production_input(self, path: Path) -> Path:
        """Reject frozen, historical, and deprecated roots as compute inputs."""

        resolved = path.resolve(strict=False)
        forbidden = (
            self.validation_root.resolve(strict=False),
            self.scrna_validation_root.resolve(strict=False),
            self.forbidden_reference_root.resolve(strict=False),
            (self.project_root / "hcc_st").resolve(strict=False),
            (self.project_root / "dat").resolve(strict=False),
            (self.project_root / "combined").resolve(strict=False),
            (self.project_root / "plan").resolve(strict=False),
            (
                self.result_root / "invalidated_frozen_input_run_20260807"
            ).resolve(strict=False),
        )
        for root in forbidden:
            if resolved == root or root in resolved.parents:
                raise ValueError(f"Forbidden production input: {resolved}")
        return resolved

    def assert_legacy_code(self, path: Path) -> Path:
        resolved = path.resolve(strict=False)
        root = self.legacy_code_root.resolve(strict=False)
        if root not in resolved.parents:
            raise ValueError(f"Legacy source must be fixed under {root}: {resolved}")
        if resolved.suffix not in {".R", ".py"}:
            raise ValueError(f"Legacy source is not executable code: {resolved}")
        return resolved

    def assert_validation_target(self, path: Path) -> Path:
        resolved = path.resolve(strict=False)
        root = self.validation_root.resolve(strict=False)
        if resolved != root and root not in resolved.parents:
            raise ValueError(f"Validation target is outside {root}: {resolved}")
        return resolved


def load_config() -> SpatialConfig:
    project_root = Path(__file__).resolve().parents[3]
    seed_text = os.environ.get("LIUF_MASTER_SEED", "20260807")
    threads_text = os.environ.get("LIUF_THREADS", "1")
    run_variant = os.environ.get(
        "LIUF_SPATIAL_RUN_VARIANT",
        "article_equivalent",
    )
    allowed_variants = {"article_equivalent", "raw_rctd_diagnostic"}
    if run_variant not in allowed_variants:
        raise ValueError(
            "LIUF_SPATIAL_RUN_VARIANT must be one of "
            f"{sorted(allowed_variants)}, got {run_variant!r}"
        )
    if re.fullmatch(r"(?:0|[1-9][0-9]*)", seed_text) is None:
        raise ValueError(
            "LIUF_MASTER_SEED must use canonical non-negative decimal "
            f"notation, got {seed_text!r}"
        )
    try:
        master_seed = int(seed_text)
    except ValueError as error:
        raise ValueError(
            f"LIUF_MASTER_SEED must be an integer, got {seed_text!r}"
        ) from error
    try:
        threads = int(threads_text)
    except ValueError as error:
        raise ValueError(
            f"LIUF_THREADS must be a positive integer, got {threads_text!r}"
        ) from error
    if not 0 <= master_seed <= R_SEED_MAX:
        raise ValueError(
            f"LIUF_MASTER_SEED must be in [0, {R_SEED_MAX}]."
        )
    if threads < 1:
        raise ValueError("LIUF_THREADS must be a positive integer.")
    config = SpatialConfig(
        project_root=project_root,
        run_variant=run_variant,
        master_seed=master_seed,
        threads=threads,
    )
    config.ensure_output_dirs()
    return config


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
