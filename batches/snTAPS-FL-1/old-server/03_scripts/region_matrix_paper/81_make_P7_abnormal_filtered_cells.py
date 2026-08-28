#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path(".")
IN = BASE / "08_qc_summary/final_single_sample/cell_filter_P7_promoter_coverage/snTAPS-FL-1.P7_with_promoter_coverage_qc.tsv"
OUT = BASE / "08_qc_summary/final_single_sample/cell_filter_P7_remove_abnormal"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(IN, sep="\t")
df["CB"] = df["CB"].astype(str)

num_cols = [
    "cluster",
    "reads_used",
    "covered_calls",
    "covered_bins",
    "modified_rate",
    "promoter_n_regions",
    "promoter_covered_calls",
    "promoter_modified_calls",
    "promoter_modified_rate",
]
for c in num_cols:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

# 全体P7的极端高覆盖阈值：只用于去掉很明显的高覆盖outlier
hi_reads = df["reads_used"].quantile(0.99)
hi_calls = df["covered_calls"].quantile(0.99)
hi_bins = df["covered_bins"].quantile(0.99)
hi_prom_calls = df["promoter_covered_calls"].quantile(0.99)

# 版本1：温和信息量过滤
# 逻辑：去掉低reads/低全局覆盖/低promoter覆盖，同时去掉极端高覆盖outlier
mild = (
    (df["reads_used"] >= 200000) &
    (df["covered_calls"] >= 400000) &
    (df["covered_bins"] >= 11000) &
    (df["promoter_n_regions"] >= 700) &
    (df["promoter_covered_calls"] >= 35000) &
    (df["reads_used"] <= hi_reads) &
    (df["covered_calls"] <= hi_calls) &
    (df["covered_bins"] <= hi_bins) &
    (df["promoter_covered_calls"] <= hi_prom_calls)
)

# 版本2：严格信息量过滤
# 逻辑：更接近cluster2–7的中位水平
strict = (
    (df["reads_used"] >= 250000) &
    (df["covered_calls"] >= 500000) &
    (df["covered_bins"] >= 12000) &
    (df["promoter_n_regions"] >= 800) &
    (df["promoter_covered_calls"] >= 40000) &
    (df["reads_used"] <= hi_reads) &
    (df["covered_calls"] <= hi_calls) &
    (df["covered_bins"] <= hi_bins) &
    (df["promoter_covered_calls"] <= hi_prom_calls)
)

# 版本3：更严格，仅用于测试是否能突出小岛
very_strict = (
    (df["reads_used"] >= 280000) &
    (df["covered_calls"] >= 550000) &
    (df["covered_bins"] >= 12500) &
    (df["promoter_n_regions"] >= 850) &
    (df["promoter_covered_calls"] >= 45000) &
    (df["reads_used"] <= hi_reads) &
    (df["covered_calls"] <= hi_calls) &
    (df["covered_bins"] <= hi_bins) &
    (df["promoter_covered_calls"] <= hi_prom_calls)
)

# 参考集：不是正式QC过滤，只是帮助确认小岛是不是更像真实细胞
island = ~df["cluster"].isin([0, 1])

filters = {
    "P7_info_mild": mild,
    "P7_info_strict": strict,
    "P7_info_very_strict": very_strict,
    "P7_island_reference_cluster2to7": island,
    "P7_info_strict_and_island": strict & island,
}

rows = []
for name, mask in filters.items():
    sub = df[mask].copy()
    sub.to_csv(OUT / f"snTAPS-FL-1.{name}.global_per_cb.tsv", sep="\t", index=False)

    row = {
        "filter": name,
        "n_cells": len(sub),
        "cluster0": int((sub["cluster"] == 0).sum()) if len(sub) else 0,
        "cluster1": int((sub["cluster"] == 1).sum()) if len(sub) else 0,
        "cluster2to7": int((~sub["cluster"].isin([0,1])).sum()) if len(sub) else 0,
        "cluster0_1_frac": float(sub["cluster"].isin([0,1]).mean()) if len(sub) else np.nan,
        "reads_used_median": sub["reads_used"].median(),
        "covered_calls_median": sub["covered_calls"].median(),
        "covered_bins_median": sub["covered_bins"].median(),
        "promoter_n_regions_median": sub["promoter_n_regions"].median(),
        "promoter_covered_calls_median": sub["promoter_covered_calls"].median(),
        "modified_rate_median": sub["modified_rate"].median(),
        "promoter_modified_rate_median": sub["promoter_modified_rate"].median(),
    }
    rows.append(row)

summary = pd.DataFrame(rows)
summary.to_csv(OUT / "P7_remove_abnormal_filter_summary.tsv", sep="\t", index=False)

print("[INFO] high outlier thresholds")
print("reads_used_q99 =", hi_reads)
print("covered_calls_q99 =", hi_calls)
print("covered_bins_q99 =", hi_bins)
print("promoter_covered_calls_q99 =", hi_prom_calls)
print()
print(summary.to_string(index=False))
