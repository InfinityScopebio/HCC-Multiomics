#!/usr/bin/env python3
"""Build patient-balanced, LR-aware pairwise CCI from production COMMOT caches."""

from __future__ import annotations

import argparse
import importlib.util
import os
import random
import sys
from pathlib import Path
from types import ModuleType

import numpy as np

from config import SpatialConfig, load_config


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def configure_runtime(config: SpatialConfig) -> None:
    os.environ["PYTHONHASHSEED"] = str(config.master_seed)
    os.environ["OMP_NUM_THREADS"] = str(config.threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(config.threads)
    os.environ["MKL_NUM_THREADS"] = str(config.threads)
    os.environ["BLIS_NUM_THREADS"] = str(config.threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(config.threads)
    os.environ["RCPP_PARALLEL_NUM_THREADS"] = str(config.threads)
    random.seed(config.master_seed)
    np.random.seed(config.master_seed)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root")
    arguments = parser.parse_args()
    config = load_config()
    configure_runtime(config)
    output_root = (
        Path(arguments.output_root).resolve()
        if arguments.output_root
        else config.cci_root
    )
    expected_output_root = config.result_root.resolve()
    if (
        output_root != expected_output_root
        and expected_output_root not in output_root.parents
    ):
        raise ValueError(f"Pairwise output must remain under {expected_output_root}")
    pairwise_path = config.assert_legacy_code(
        config.legacy_code_root / "cci_lr_aware_pairwise.py"
    )
    pairwise = load_module("cci_lr_aware_pairwise", pairwise_path)
    pairwise.EXPECTED_CACHE_VERSION = config.active_commot_cache_version
    os.chdir(config.project_root)

    display = {
        "HSC": "HSC",
        "Pericyte": "Pericyte",
        "Treg": "Treg",
        "IL1B_Mp": "IL1B+Mφ",
        "SPP1_Mp": "SPP1+Mφ",
        "TREM2_Mp": "TREM2+Mφ",
        "CXCL10_Mp": "CXCL10+Mφ",
        "Hep_C4": "Hep_C4",
    }
    pairwise_config = pairwise.PairwiseConfig(
        commot_dir=config.commot_root,
        output_dir=output_root,
        sections=config.tumor_sections,
        cm02_types=config.cm02_types,
        display_labels=display,
        mean_profile_file=config.cci_root / "rctd_reference_lr_mean_cpm.csv",
        pct_profile_file=config.cci_root / "rctd_reference_lr_pct.csv",
        min_detection_fraction=0.05,
        top_pathways_per_edge=5,
    )
    results = pairwise.run_pairwise_analysis(pairwise_config)
    pairwise.save_pairwise_results(results, pairwise_config)
    print(
        f"Pairwise CCI: {len(results.lr_edges):,} LR rows; "
        f"{len(results.edges):,} directed cell-state edges."
    )


if __name__ == "__main__":
    main()
