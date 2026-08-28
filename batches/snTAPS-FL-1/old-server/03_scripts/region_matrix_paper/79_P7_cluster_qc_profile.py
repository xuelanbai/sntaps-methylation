#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path(".")
IN = BASE / "08_qc_summary/final_single_sample/cell_filter_P7_promoter_coverage/snTAPS-FL-1.P7_with_promoter_coverage_qc.tsv"
OUT = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_cluster_qc_profile"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(IN, sep="\t")
df["CB"] = df["CB"].astype(str)
df["cluster"] = pd.to_numeric(df["cluster"], errors="coerce").astype("Int64")

metrics = [
    "reads_used",
    "covered_calls",
    "covered_bins",
    "modified_rate",
    "promoter_n_regions",
    "promoter_covered_calls",
    "promoter_modified_calls",
    "promoter_modified_rate",
]

metrics = [c for c in metrics if c in df.columns]

for c in metrics:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# 全体P7阈值，用来判断每个cluster是否富集极端细胞
thresholds = {}
for c in metrics:
    thresholds[f"{c}_q01"] = df[c].quantile(0.01)
    thresholds[f"{c}_q05"] = df[c].quantile(0.05)
    thresholds[f"{c}_q10"] = df[c].quantile(0.10)
    thresholds[f"{c}_q50"] = df[c].quantile(0.50)
    thresholds[f"{c}_q90"] = df[c].quantile(0.90)
    thresholds[f"{c}_q95"] = df[c].quantile(0.95)
    thresholds[f"{c}_q99"] = df[c].quantile(0.99)

pd.Series(thresholds).to_csv(OUT / "P7_global_metric_thresholds.tsv", sep="\t", header=["value"])

rows = []
for cl, sub in df.groupby("cluster", dropna=False):
    row = {
        "cluster": cl,
        "n_cells": len(sub),
        "frac_cells": len(sub) / len(df),
    }

    for c in metrics:
        row[f"{c}_mean"] = sub[c].mean()
        row[f"{c}_median"] = sub[c].median()
        row[f"{c}_q10"] = sub[c].quantile(0.10)
        row[f"{c}_q90"] = sub[c].quantile(0.90)
        row[f"{c}_q95"] = sub[c].quantile(0.95)

        q05 = thresholds[f"{c}_q05"]
        q95 = thresholds[f"{c}_q95"]
        q99 = thresholds[f"{c}_q99"]

        row[f"{c}_low5_frac"] = (sub[c] <= q05).mean()
        row[f"{c}_high95_frac"] = (sub[c] >= q95).mean()
        row[f"{c}_high99_frac"] = (sub[c] >= q99).mean()

    rows.append(row)

summary = pd.DataFrame(rows).sort_values("cluster")
summary.to_csv(OUT / "P7_cluster_qc_summary.tsv", sep="\t", index=False)

# 做一个更易读的median矩阵
median_mat = summary[["cluster", "n_cells", "frac_cells"] + [f"{c}_median" for c in metrics]]
median_mat.to_csv(OUT / "P7_cluster_qc_median_table.tsv", sep="\t", index=False)

# z-score化后的cluster median，用于快速看哪个cluster偏高/偏低
z_rows = []
for c in metrics:
    vals = summary[f"{c}_median"].astype(float)
    z = (vals - vals.mean()) / vals.std(ddof=0) if vals.std(ddof=0) > 0 else vals * 0
    for cl, zv in zip(summary["cluster"], z):
        z_rows.append({"metric": c, "cluster": cl, "cluster_median_zscore": zv})

zdf = pd.DataFrame(z_rows)
zdf.to_csv(OUT / "P7_cluster_qc_median_zscore.tsv", sep="\t", index=False)

print("[INFO] input:", IN)
print("[INFO] output:", OUT)
print("[INFO] n_cells:", len(df))
print("[INFO] metrics:", metrics)
print(summary[["cluster", "n_cells", "frac_cells"] + [f"{c}_median" for c in metrics]].to_string(index=False))
