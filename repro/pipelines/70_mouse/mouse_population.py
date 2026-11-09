"""Read the verified animal population before any endpoint calculations."""

from __future__ import annotations

import csv
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import openpyxl


@dataclass(frozen=True)
class PopulationConfig:
    workbook: Path
    exclusion_csv: Path
    mirrored_exclusion_csv: Path
    groups: Tuple[str, ...] = tuple("ABCDEFG")
    n_per_arm: int = 5
    expected_external_records: Tuple[Tuple[str, int], ...] = (("E", 21),)
    confirmation_date: str = "2026-09-10"


@dataclass(frozen=True)
class Animal:
    group: str
    sample_id: int
    terminal_body_weight_g: float
    tumour_weight_g: float
    worksheet_row: int


@dataclass(frozen=True)
class Population:
    animals: Tuple[Animal, ...]
    external_records: Tuple[Tuple[str, int, int], ...]
    exclusion_records: Tuple[Dict[str, str], ...]
    worksheet_record_count: int


def default_config(project_root: Path) -> PopulationConfig:
    """Use the immutable raw workbook and the two versioned population ledgers."""
    return PopulationConfig(
        workbook=project_root / "repo/data/private/raw/mouse"
        / "小鼠样本图片及送样单/终末取材数据.xlsx",
        exclusion_csv=project_root / "repo/data/reference/mouse"
        / "Fig7_analysis_population_exclusions.csv",
        mirrored_exclusion_csv=project_root / "plan/main_figure_final/Fig7"
        / "Fig7_analysis_population_exclusions.csv",
    )


def read_exclusion_records(config: PopulationConfig) -> Tuple[Dict[str, str], ...]:
    """Reject drift between the two ledgers or unsupported population changes."""
    ledgers: List[Tuple[Dict[str, str], ...]] = []
    required = {"group", "sample_id", "analysis_status", "reason", "basis"}
    for path in (config.exclusion_csv, config.mirrored_exclusion_csv):
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(f"Invalid animal exclusion ledger: {path}")
            ledgers.append(tuple(dict(row) for row in reader))
    if ledgers[0] != ledgers[1]:
        raise ValueError("The two animal analysis-population ledgers differ.")
    records = ledgers[0]
    keys = [(row["group"], int(row["sample_id"])) for row in records]
    if len(keys) != len(set(keys)) or set(keys) != set(config.expected_external_records):
        raise ValueError("The animal ledger contains unexpected external records.")
    for row in records:
        if row["analysis_status"] != "exclude":
            raise ValueError("External records must be excluded from all analyses.")
        if "not a sample from this experiment" not in row["reason"].lower():
            raise ValueError("The external-record reason is missing or outdated.")
        if config.confirmation_date not in row["basis"]:
            raise ValueError("The user confirmation date is missing.")
    return records


def load_population(config: PopulationConfig) -> Population:
    """Exclude external records before reading endpoint values or calculating means."""
    ledger = read_exclusion_records(config)
    external_keys = {(row["group"], int(row["sample_id"])) for row in ledger}
    workbook = openpyxl.load_workbook(
        config.workbook, data_only=True, read_only=True
    )
    animals: List[Animal] = []
    external: List[Tuple[str, int, int]] = []
    raw_keys: List[Tuple[str, int]] = []
    current_group = ""
    try:
        sheet = workbook.worksheets[0]
        header = next(sheet.iter_rows(values_only=True))
        if tuple(header[index] for index in (1, 3, 4, 5)) != (
            "分组", "样本编号", "小鼠体重(g)", "肝部肿瘤（g）"
        ):
            raise ValueError("Unexpected terminal-collection workbook columns.")
        for row_number, row in enumerate(
            sheet.iter_rows(min_row=2, values_only=True), start=2
        ):
            if row[1]:
                current_group = str(row[1]).replace("组", "").strip()
            if row[3] is None:
                continue
            if current_group not in config.groups:
                raise ValueError(f"Unknown group at worksheet row {row_number}.")
            sample_id = int(row[3])
            if float(row[3]) != sample_id:
                raise ValueError(f"Non-integer animal ID at row {row_number}.")
            key = (current_group, sample_id)
            raw_keys.append(key)
            if key in external_keys:
                external.append((current_group, sample_id, row_number))
                continue
            body, tumour = float(row[4]), float(row[5])
            if not all(math.isfinite(value) and value > 0 for value in (body, tumour)):
                raise ValueError(f"Invalid terminal endpoint in retained animal {key}.")
            animals.append(Animal(current_group, sample_id, body, tumour, row_number))
    finally:
        workbook.close()

    if len(raw_keys) != len(set(raw_keys)) or len(raw_keys) != len(
        {sample_id for _, sample_id in raw_keys}
    ):
        raise ValueError("Animal IDs must be unique in the source workbook.")
    if {(group, sample_id) for group, sample_id, _ in external} != external_keys:
        raise ValueError("Each confirmed external record must occur exactly once.")
    counts = Counter(animal.group for animal in animals)
    if counts != Counter({group: config.n_per_arm for group in config.groups}):
        raise ValueError(f"Expected five experimental animals per arm, found {counts}.")
    if len(raw_keys) != len(animals) + len(external_keys):
        raise ValueError("Source-to-analysis animal accounting failed.")
    return Population(tuple(animals), tuple(external), ledger, len(raw_keys))
