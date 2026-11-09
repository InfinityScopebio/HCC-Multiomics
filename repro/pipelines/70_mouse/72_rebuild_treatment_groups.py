#!/usr/bin/env python3
"""Regenerate Fig7_treatment_groups.csv from the raw terminal-collection workbook.

Analysis set: n = 5 per arm, 35 mice. The versioned analysis-population record
excludes external record E-21 (not from this experiment; confirmed 2026-09-10).
The shared reader ensures that the summary
table, the displayed panel and the pairwise audit table describe the same mice.
"""
import csv
import os
import random
from dataclasses import replace
from pathlib import Path

import numpy as np

from scipy import stats

from mouse_population import default_config, load_population

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RESULT_DIR = os.environ.get(
    "LIUF_MOUSE_RESULT_DIR", os.path.join(ROOT, "repo", "results", "70_mouse")
)
os.makedirs(RESULT_DIR, exist_ok=True)
MASTER_SEED = int(os.environ.get("LIUF_MASTER_SEED", "20260807"))
random.seed(MASTER_SEED)
np.random.seed(MASTER_SEED)
XLSX = os.path.join(
    ROOT,
    "repo", "data", "private", "raw", "mouse",
    "小鼠样本图片及送样单", "终末取材数据.xlsx",
)
OUT = os.path.join(RESULT_DIR, "Fig7_treatment_groups.csv")

G = list("ABCDEFG")
EN = {"A": "Vehicle", "B": "Hydroxynidone", "C": "Lenvatinib", "D": "Trametinib",
      "E": "Lenvatinib + Trametinib", "F": "Lenvatinib + Hydroxynidone",
      "G": "Lenvatinib + Trametinib + Hydroxynidone"}
CN = {"A": "对照组（溶媒）", "B": "羟尼酮组", "C": "仑伐替尼组", "D": "曲美替尼组",
      "E": "仑伐替尼 + 曲美替尼", "F": "仑伐替尼 + 羟尼酮", "G": "三药联用"}
DOSE = {"A": "0.5% CMC-Na", "B": "50 mg/kg (50% of 100)", "C": "8 mg/kg",
        "D": "0.5 mg/kg (50% of 1)", "E": "8 + 0.5 mg/kg", "F": "8 + 50 mg/kg",
        "G": "8 + 0.5 + 50 mg/kg"}

# All endpoints share the verified 35-animal population and dated ledgers.
population_config = replace(default_config(Path(ROOT)), workbook=Path(XLSX))
population = load_population(population_config)
rec = {
    group: [
        (animal.sample_id, animal.terminal_body_weight_g, animal.tumour_weight_g)
        for animal in population.animals if animal.group == group
    ]
    for group in G
}
n_total = len(population.animals)

# ---------------- summarise ----------------
veh = np.mean([t for _, _, t in rec["A"]])
rows = []
for k in G:
    ids = [s for s, _, _ in rec[k]]
    body = np.array([b for _, b, _ in rec[k]])
    tum = np.array([t for _, _, t in rec[k]])
    rows.append({
        "Group": k, "Treatment_EN": EN[k], "Treatment_CN": CN[k], "Dose": DOSE[k],
        "n": len(ids), "Sample_IDs": " ".join(str(i) for i in ids),
        "Tumour_g_mean": round(float(tum.mean()), 2),
        "Tumour_g_SD": round(float(tum.std(ddof=1)), 2),
        "Tumour_range": f"{tum.min():.1f}–{tum.max():.1f}",
        "Body_g_mean": round(float(body.mean()), 1),
        "vs_Vehicle": f"{round(100 * tum.mean() / veh):d}%",
    })

with open(OUT, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

F, p = stats.f_oneway(*[np.array([t for _, _, t in rec[k]]) for k in G])
print(f"wrote {OUT}  (n = {n_total}, 5 per arm)")
print(f"one-way ANOVA F({len(G)-1},{n_total-len(G)}) = {F:.2f}  P = {p:.3g}")
for r in rows:
    print(f"  {r['Group']} {r['Treatment_EN']:<40s} n={r['n']} "
          f"mean={r['Tumour_g_mean']:.2f} SD={r['Tumour_g_SD']:.2f} "
          f"range={r['Tumour_range']} vs_veh={r['vs_Vehicle']}")
