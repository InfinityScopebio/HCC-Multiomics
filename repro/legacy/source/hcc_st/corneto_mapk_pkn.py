"""Curated MAPK prior-network definition shared by the CORNETO runner."""

from __future__ import annotations

import re


def split_genes(value: str) -> set[str]:
    """Split CellChat complex or overlap-gene fields into symbols."""
    return {
        token.strip().upper()
        for token in re.split(r"[_+&|:;]", str(value))
        if token.strip()
    }


def build_mapk_pkn() -> list[tuple[str, int, str]]:
    """Return the literature-curated receptor-to-ERK prior network."""
    edges = [
        ("FGFR2", 1, "FRS2"),
        ("FGFR2", 1, "SHC1"),
        ("TGFBR2", 1, "TGFBR1"),
        ("TGFBR1", 1, "SHC1"),
        ("FRS2", 1, "GRB2"),
        ("SHC1", 1, "GRB2"),
        ("GRB2", 1, "SOS1"),
        ("SOS1", 1, "HRAS"),
        ("SOS1", 1, "NRAS"),
        ("SOS1", 1, "KRAS"),
        ("HRAS", 1, "RAF1"),
        ("HRAS", 1, "BRAF"),
        ("NRAS", 1, "RAF1"),
        ("NRAS", 1, "BRAF"),
        ("KRAS", 1, "RAF1"),
        ("KRAS", 1, "BRAF"),
        ("KRAS", 1, "ARAF"),
        ("RAF1", 1, "MAP2K1"),
        ("RAF1", 1, "MAP2K2"),
        ("BRAF", 1, "MAP2K1"),
        ("BRAF", 1, "MAP2K2"),
        ("ARAF", 1, "MAP2K1"),
        ("MAP2K1", 1, "MAPK1"),
        ("MAP2K1", 1, "MAPK3"),
        ("MAP2K2", 1, "MAPK1"),
        ("MAP2K2", 1, "MAPK3"),
        ("MAPK1", 1, "MYC"),
        ("MAPK1", 1, "ELK1"),
        ("MAPK1", 1, "FOS"),
        ("MAPK1", 1, "JUN"),
        ("MAPK1", 1, "ETS1"),
        ("MAPK1", 1, "CREB1"),
        ("MAPK3", 1, "MYC"),
        ("MAPK3", 1, "ELK1"),
        ("MAPK3", 1, "FOS"),
        ("MAPK3", 1, "ETS1"),
    ]
    return list(dict.fromkeys(edges))
