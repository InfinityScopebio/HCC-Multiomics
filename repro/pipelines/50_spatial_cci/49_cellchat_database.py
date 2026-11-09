#!/usr/bin/env python3
"""Materialize the COMMOT-bundled CellChat database with provenance."""

from __future__ import annotations

import json
import os
import random
from pathlib import Path

import commot
import numpy as np
import pandas as pd

from config import load_config, sha256_file


def main() -> None:
    config = load_config()
    os.environ["PYTHONHASHSEED"] = str(config.master_seed)
    random.seed(config.master_seed)
    np.random.seed(config.master_seed)
    package_root = Path(commot.__file__).resolve().parent
    database_path = (
        package_root
        / "preprocessing/_data/LRdatabase/CellChat/CellChatDB.ligrec.human.csv"
    )
    config.assert_production_input(database_path)
    if not database_path.is_file():
        raise FileNotFoundError(database_path)
    database = pd.read_csv(database_path, index_col=0)
    database.columns = ["ligand", "receptor", "pathway", "annotation"]
    output_path = config.commot_root / "cellchat_ligrec_full.csv"
    database.to_csv(output_path, index=False)
    provenance = {
        "source": str(database_path),
        "source_sha256": sha256_file(database_path),
        "output": str(output_path),
        "output_sha256": sha256_file(output_path),
        "commot_version": getattr(commot, "__version__", "unknown"),
        "rows": len(database),
        "usage": "production package database",
    }
    (config.commot_root / "cellchat_database_provenance.json").write_text(
        json.dumps(provenance, indent=2),
        encoding="utf-8",
    )
    print(f"CellChat LR database: {len(database):,} rows.")


if __name__ == "__main__":
    main()
