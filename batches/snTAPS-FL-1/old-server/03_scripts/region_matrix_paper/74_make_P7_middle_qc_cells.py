#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path(".")
P7 = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"
OUT = BASE / "08_qc_summary/final_single_sample/cell_filter_P7_middle_qc"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(P7, sep="\t")
df["CB"] = df["CB"].astype(str)

qc_cols = ["reads_used", "covered_calls", "covered_bins", "modified_rate"]
for c in qc_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

print("[INFO] P7 cells:", len(df), flush=True)

# A: 去低reads一半，同时去掉最高5%覆盖异常
d = df.copy()
reads_median = d["reads_used"].median()
m = d["reads_used"] >= reads_median
for c in ["reads_used", "covered_calls", "covered_bins"]:
    hi = d[c].quantile(0.95)
    m &= d[c] <= hi
# modified_rate只去极端两端
m &= d["modified_rate"].between(d["modified_rate"].quantile(0.05), d["modified_rate"].quantile(0.95))
A = d[m].copy()
A.to_csv(OUT / "snTAPS-FL-1.P7_keep_reads_ge_median_dropHi95.global_per_cb.tsv", sep="\t", index=False)

# B: 所有核心QC在10%-90%
m = pd.Series(True, index=df.index)
for c in qc_cols:
    m &= df[c].between(df[c].quantile(0.10), df[c].quantile(0.90))
B = df[m].copy()
B.to_csv(OUT / "snTAPS-FL-1.P7_middle80_qc.global_per_cb.tsv", sep="\t", index=False)

# C: 所有核心QC在25%-75%
m = pd.Series(True, index=df.index)
for c in qc_cols:
    m &= df[c].between(df[c].quantile(0.25), df[c].quantile(0.75))
C = df[m].copy()
C.to_csv(OUT / "snTAPS-FL-1.P7_middle50_qc.global_per_cb.tsv", sep="\t", index=False)

summary = []
for name, sub in [
    ("P7_keep_reads_ge_median_dropHi95", A),
    ("P7_middle80_qc", B),
    ("P7_middle50_qc", C),
]:
    row = {"filter": name, "n_cells": len(sub)}
    if "cluster" in sub.columns:
        row["cluster0"] = int((sub["cluster"] == 0).sum())
        row["cluster1"] = int((sub["cluster"] == 1).sum())
        row["cluster0_1"] = row["cluster0"] + row["cluster1"]
        row["cluster0_1_frac"] = row["cluster0_1"] / len(sub) if len(sub) else np.nan
    for c in qc_cols:
        row[f"{c}_median"] = sub[c].median()
    summary.append(row)

pd.DataFrame(summary).to_csv(OUT / "P7_middle_qc_filter_summary.tsv", sep="\t", index=False)

print("[INFO] wrote:", OUT, flush=True)
print(pd.DataFrame(summary), flush=True)
