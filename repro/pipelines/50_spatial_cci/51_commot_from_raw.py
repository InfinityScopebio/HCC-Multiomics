#!/usr/bin/env python3
"""Run COMMOT on raw Visium expression at production maxCM == CM02 spots."""

from __future__ import annotations

import argparse
import json
import os
import random
from pathlib import Path
from typing import Any

import commot as ct
import numpy as np
import pandas as pd
import pyreadr
import scanpy as sc
from scipy.spatial.distance import pdist

from config import SpatialConfig, derive_r_seed, load_config, sha256_file


HETEROMERIC_RULE = "min"
HETEROMERIC_DELIMITER = "_"
WEIGHT_METADATA = {"barcode", "row", "col", "imgrow", "imgcol"}


def section_seed(config: SpatialConfig, section: str) -> int:
    """Derive a stable seed independent of requested section order."""

    return derive_r_seed(
        config.master_seed,
        config.all_sections.index(section) + 1,
    )


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


def read_positions(section_dir: Path) -> pd.DataFrame:
    path = section_dir / "spatial/tissue_positions_list.csv"
    positions = pd.read_csv(path, header=None)
    positions.columns = (
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row",
        "pxl_col",
    )
    return positions.set_index("barcode")


def load_raw_visium(config: SpatialConfig, section: str) -> sc.AnnData:
    section_dir = config.assert_production_input(config.visium_dir(section))
    matrix_path = section_dir / "filtered_feature_bc_matrix.h5"
    adata = sc.read_10x_h5(matrix_path)
    adata.var_names_make_unique()
    positions = read_positions(section_dir)
    common = [
        barcode
        for barcode in adata.obs_names
        if barcode in positions.index and positions.at[barcode, "in_tissue"] == 1
    ]
    adata = adata[common].copy()
    selected_positions = positions.loc[adata.obs_names]
    adata.obs["in_tissue"] = selected_positions["in_tissue"].astype(int)
    adata.obsm["spatial"] = selected_positions[
        ["pxl_col", "pxl_row"]
    ].to_numpy(dtype=float)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def load_weights(config: SpatialConfig, section: str) -> pd.DataFrame:
    path = config.assert_production_input(config.active_weight_rds(section))
    objects = pyreadr.read_r(str(path))
    if not objects:
        raise ValueError(f"{section}: active RCTD RDS contains no data frame")
    weights = next(iter(objects.values()))
    if "barcode" not in weights.columns:
        raise ValueError(f"{section}: historical RCTD weights lack barcode")
    weights = weights.set_index("barcode")
    if weights.index.has_duplicates:
        raise ValueError(f"{section}: duplicate RCTD barcodes")
    return weights


def load_recovery(config: SpatialConfig) -> pd.DataFrame:
    path = config.assert_production_input(
        config.cm_recovery_with_barcode_rds
    )
    objects = pyreadr.read_r(str(path))
    if not objects:
        raise ValueError("CM recovery RDS contains no data frame")
    recovery = next(iter(objects.values()))
    required = {"section", "barcode", "maxCM"}
    if not required.issubset(recovery.columns):
        raise ValueError(f"CM recovery lacks {sorted(required)}")
    return recovery


def cellchat_database(config: SpatialConfig) -> pd.DataFrame:
    path = config.assert_production_input(
        config.commot_root / "cellchat_ligrec_full.csv"
    )
    database = pd.read_csv(path)
    return database[["ligand", "receptor", "pathway", "annotation"]]


def numeric_fingerprint(
    config: SpatialConfig,
    adata: sc.AnnData,
    section: str,
    *,
    cache_path: Path | None = None,
) -> dict[str, Any]:
    keys = sorted(
        key
        for key in adata.obsp
        if key.startswith("commot-CellChat-")
    )
    matrix_sum = 0.0
    matrix_square_sum = 0.0
    matrix_nnz = 0
    for key in keys:
        matrix = adata.obsp[key]
        values = matrix.data if hasattr(matrix, "data") else np.asarray(matrix).ravel()
        matrix_sum += float(np.asarray(values, dtype=float).sum())
        matrix_square_sum += float(np.square(np.asarray(values, dtype=float)).sum())
        matrix_nnz += int(np.count_nonzero(values))
    return {
        "section": section,
        "n_spots": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_commot_matrices": len(keys),
        "commot_nnz": matrix_nnz,
        "commot_sum": matrix_sum,
        "commot_square_sum": matrix_square_sum,
        "rctd_sum": float(np.asarray(adata.obsm["rctd_all"]).sum()),
        "seed": int(adata.uns.get("section_seed", section_seed(config, section))),
        "threads": int(adata.uns.get("threads", config.threads)),
        "filtered_lr_pairs": int(adata.uns.get("commot_filtered_lr_count", -1)),
        "cache_sha256": (
            sha256_file(cache_path)
            if cache_path is not None and cache_path.is_file()
            else ""
        ),
    }


def run_section(
    config: SpatialConfig,
    section: str,
    output_root: Path,
    database: pd.DataFrame,
) -> tuple[Path, dict[str, Any]]:
    current_seed = section_seed(config, section)
    random.seed(current_seed)
    np.random.seed(current_seed)
    adata = load_raw_visium(config, section)
    weights = load_weights(config, section)
    recovery = load_recovery(config)
    common = [barcode for barcode in adata.obs_names if barcode in weights.index]
    adata = adata[common].copy()
    weights = weights.loc[common].copy()
    selected = recovery[
        (recovery["section"] == section) & (recovery["maxCM"] == "CM02")
    ]["barcode"].astype(str)
    selected_set = set(selected)
    selected_barcodes = [
        barcode for barcode in adata.obs_names if barcode in selected_set
    ]
    if len(selected_barcodes) != len(selected_set):
        missing = sorted(selected_set.difference(selected_barcodes))
        raise ValueError(
            f"{section}: {len(missing)} production CM02 barcodes missing from raw Visium"
        )
    adata = adata[selected_barcodes].copy()
    weights = weights.loc[selected_barcodes]
    if adata.n_obs < 30:
        raise ValueError(f"{section}: fewer than 30 maxCM == CM02 spots")

    coordinates = adata.obsm["spatial"].astype(float)
    distances = pdist(coordinates)
    positive_distances = distances[distances > 0]
    if len(positive_distances) == 0:
        raise ValueError(f"{section}: no positive spatial distance")
    adata.obsm["spatial"] = coordinates * (
        100.0 / float(positive_distances.min())
    )

    weight_columns = [
        column for column in weights.columns if column not in WEIGHT_METADATA
    ]
    missing_cm02 = set(config.cm02_types).difference(weight_columns)
    if missing_cm02:
        raise ValueError(
            f"{section}: RCTD weights lack {sorted(missing_cm02)}"
        )
    adata.obsm["rctd_all"] = weights[weight_columns].to_numpy(dtype=float)
    adata.uns["rctd_all_cols"] = weight_columns
    adata.obsm["rctd_cm02"] = weights[list(config.cm02_types)].to_numpy(
        dtype=float
    )
    adata.uns["rctd_cm02_cols"] = list(config.cm02_types)
    adata.uns["cm02_selection_source"] = (
        "production_CM_recover_barcode_join_maxCM_equals_CM02"
    )
    adata.uns["cci_cache_version"] = config.active_commot_cache_version
    adata.uns["rctd_source"] = config.rctd_source
    adata.uns["article_equivalent"] = config.article_equivalent
    adata.uns["commot_heteromeric"] = True
    adata.uns["commot_heteromeric_rule"] = HETEROMERIC_RULE
    adata.uns["commot_heteromeric_delimiter"] = HETEROMERIC_DELIMITER
    adata.uns["master_seed"] = config.master_seed
    adata.uns["section_seed"] = current_seed
    adata.uns["threads"] = config.threads

    filtered_database = ct.pp.filter_lr_database(
        database,
        adata,
        heteromeric=True,
        heteromeric_delimiter=HETEROMERIC_DELIMITER,
        heteromeric_rule=HETEROMERIC_RULE,
        min_cell_pct=0.01,
    )
    if filtered_database.empty:
        raise RuntimeError(f"{section}: no CellChat LR pairs after filtering")
    adata.uns["commot_filtered_lr_count"] = int(len(filtered_database))
    ct.tl.spatial_communication(
        adata,
        database_name="CellChat",
        df_ligrec=filtered_database,
        pathway_sum=True,
        heteromeric=True,
        heteromeric_rule=HETEROMERIC_RULE,
        heteromeric_delimiter=HETEROMERIC_DELIMITER,
        dis_thr=200.0,
    )
    output_root.mkdir(parents=True, exist_ok=True)
    output_path = output_root / f"commot_{section}.h5ad"
    temporary_path = output_path.with_suffix(".h5ad.tmp")
    adata.write_h5ad(temporary_path)
    os.replace(temporary_path, output_path)
    fingerprint = numeric_fingerprint(
        config,
        adata,
        section,
        cache_path=output_path,
    )
    return output_path, fingerprint


def compare_fingerprints(
    first: dict[str, Any],
    second: dict[str, Any],
) -> dict[str, Any]:
    numeric_fields = (
        "n_spots",
        "n_genes",
        "n_commot_matrices",
        "commot_nnz",
        "commot_sum",
        "commot_square_sum",
        "rctd_sum",
        "seed",
        "threads",
        "filtered_lr_pairs",
    )
    deltas = {
        field: abs(float(first[field]) - float(second[field]))
        for field in numeric_fields
    }
    maximum = max(deltas.values())
    return {
        "step": "COMMOT",
        "section": first["section"],
        "repeat_1": json.dumps(first, sort_keys=True),
        "repeat_2": json.dumps(second, sort_keys=True),
        "maximum_numeric_delta": maximum,
        "status": "exact" if maximum == 0 else "blocked_non_deterministic",
        "seed": int(first["seed"]),
        "threads": int(first["threads"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sections", nargs="+")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--repeat-check", action="store_true")
    arguments = parser.parse_args()
    config = load_config()
    configure_runtime(config)
    sections = tuple(arguments.sections or config.tumor_sections)
    invalid = set(sections).difference(config.tumor_sections)
    if invalid:
        raise ValueError(f"COMMOT is tumor-only; invalid sections: {sorted(invalid)}")
    database = cellchat_database(config)
    fingerprints: list[dict[str, Any]] = []
    for section in sections:
        output_path = config.commot_cache(section)
        if output_path.is_file() and not arguments.force:
            cached = sc.read_h5ad(output_path)
            if (
                str(cached.uns.get("cci_cache_version", ""))
                == config.active_commot_cache_version
            ):
                fingerprints.append(
                    numeric_fingerprint(
                        config,
                        cached,
                        section,
                        cache_path=output_path,
                    )
                )
                print(f"{section}: cached production COMMOT")
                continue
        _, fingerprint = run_section(
            config,
            section,
            config.commot_root,
            database,
        )
        fingerprints.append(fingerprint)
        print(f"{section}: {fingerprint['n_spots']} CM02-positive spots")
    pd.DataFrame(fingerprints).to_csv(
        config.commot_root / "commot_numeric_fingerprint.tsv",
        sep="\t",
        index=False,
    )

    if arguments.repeat_check:
        section = sections[0]
        repeat_root = config.result_root / "repeatability/commot_repeat"
        _, first = run_section(config, section, repeat_root / "run1", database)
        _, second = run_section(config, section, repeat_root / "run2", database)
        comparison = compare_fingerprints(first, second)
        pd.DataFrame([comparison]).to_csv(
            config.audit_root / "spatial_commot_repeatability.tsv",
            sep="\t",
            index=False,
        )
        if comparison["status"] != "exact":
            raise RuntimeError("COMMOT repeatability check failed.")


if __name__ == "__main__":
    main()
