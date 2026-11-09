#!/usr/bin/env python3
"""Run the fixed LR-semantic mapping against production CCI outputs."""

from __future__ import annotations

import argparse
import importlib.util
import os
import random
import sys
from pathlib import Path

import numpy as np

from config import load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root")
    arguments = parser.parse_args()
    config = load_config()
    os.environ["PYTHONHASHSEED"] = str(config.master_seed)
    random.seed(config.master_seed)
    np.random.seed(config.master_seed)
    source = config.assert_legacy_code(
        config.legacy_code_root / "run_CCI_B10_lr_overlap_program_mapping.py"
    )
    spec = importlib.util.spec_from_file_location("fixed_lr_overlap", source)
    if spec is None or spec.loader is None:
        raise ImportError(source)
    module = importlib.util.module_from_spec(spec)
    sys.modules["fixed_lr_overlap"] = module
    spec.loader.exec_module(module)
    os.chdir(config.project_root)
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
        raise ValueError(f"River output must remain under {expected_output_root}")
    module.OUT = str(output_root)
    module.PROGRAM_GENES = str(
        config.cci_root / "target_sets_4pathway_expressed_long.csv"
    )
    module.LR_AWARE_FLOW_META = str(
        config.cci_root / "true_pairwise_lr_edges.csv"
    )
    module.main()
    print(f"LR river outputs: {config.cci_root}")


if __name__ == "__main__":
    main()
