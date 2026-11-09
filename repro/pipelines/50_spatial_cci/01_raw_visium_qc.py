#!/usr/bin/env python3
"""Inspect all raw Visium matrices without loading any derived spatial object."""

from __future__ import annotations

import json
import os
import random
from typing import Any

import h5py
import numpy as np
import pandas as pd

from config import SpatialConfig, load_config, sha256_file


def decode(values: np.ndarray) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def inspect_section(config: SpatialConfig, section: str) -> dict[str, Any]:
    section_root = config.assert_production_input(config.visium_dir(section))
    matrix_path = section_root / "filtered_feature_bc_matrix.h5"
    positions_path = section_root / "spatial/tissue_positions_list.csv"
    scale_path = section_root / "spatial/scalefactors_json.json"
    with h5py.File(matrix_path, "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(int(value) for value in matrix["shape"][...])
        barcodes = decode(matrix["barcodes"][...])
        data = matrix["data"]
        nonzero = int(data.shape[0])
        count_sum = int(np.asarray(data[...], dtype=np.int64).sum())
    positions = pd.read_csv(positions_path, header=None)
    positions.columns = (
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row",
        "pxl_col",
    )
    with scale_path.open(encoding="utf-8") as handle:
        scale = json.load(handle)
    barcode_set = set(barcodes)
    position_set = set(positions["barcode"].astype(str))
    return {
        "section": section,
        "region": section[-1],
        "is_tumor_scope": section in config.tumor_sections,
        "n_genes": shape[0],
        "n_filtered_spots": shape[1],
        "matrix_nnz": nonzero,
        "raw_count_sum": count_sum,
        "n_position_rows": len(positions),
        "n_in_tissue_position_rows": int(positions["in_tissue"].eq(1).sum()),
        "barcodes_in_positions": len(barcode_set & position_set),
        "missing_position_barcodes": len(barcode_set - position_set),
        "spot_diameter_fullres": scale.get("spot_diameter_fullres"),
        "matrix_sha256": sha256_file(matrix_path),
        "positions_sha256": sha256_file(positions_path),
        "status": (
            "exact_raw_qc"
            if shape[1] == len(barcode_set & position_set)
            else "failed_barcode_alignment"
        ),
    }


def main() -> None:
    config = load_config()
    os.environ["PYTHONHASHSEED"] = str(config.master_seed)
    random.seed(config.master_seed)
    np.random.seed(config.master_seed)
    rows = [inspect_section(config, section) for section in config.all_sections]
    table = pd.DataFrame(rows)
    result_path = config.result_root / "raw_qc/visium_section_qc.tsv"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(result_path, sep="\t", index=False)
    table.to_csv(
        config.audit_root / "spatial_raw_visium_qc.tsv",
        sep="\t",
        index=False,
    )
    if not table["status"].eq("exact_raw_qc").all():
        raise RuntimeError("At least one raw Visium section failed barcode alignment.")
    print(
        f"Raw Visium QC: {len(table)} sections, "
        f"{int(table['n_filtered_spots'].sum()):,} filtered spots."
    )


if __name__ == "__main__":
    main()
