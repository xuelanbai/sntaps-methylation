#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

IN = Path("08_qc_summary/final_single_sample/cell_filter_P7_promoter_coverage/snTAPS-FL-1.P7_with_promoter_coverage_qc.tsv")
OUT = Path("08_qc_summary/final_single_sample/cell_filter_P7_QCpass")
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

# 去掉极端高覆盖，避免doublet/multiplet-like outlier
q99_reads = df["reads_used"].quantile(0.99)
q99_calls = df["covered_calls"].quantile(0.99)
q99_bins = df["covered_bins"].quantile(0.99)
q99_prom = df["promoter_covered_calls"].quantile(0.99)

strict = (
    (df["reads_used"] >= 250000) &
    (df["covered_calls"] >= 500000) &
    (df["covered_bins"] >= 12000) &
    (df["promoter_n_regions"] >= 850) &
    (df["promoter_covered_calls"] >= 40000) &
    (df["reads_used"] <= q99_reads) &
    (df["covered_calls"] <= q99_calls) &
    (df["covered_bins"] <= q99_bins) &
    (df["promoter_covered_calls"] <= q99_prom)
)

very_strict = (
    (df["reads_used"] >= 280000) &
    (df["covered_calls"] >= 550000) &
    (df["covered_bins"] >= 12500) &
    (df["promoter_n_regions"] >= 900) &
    (df["promoter_covered_calls"] >= 45000) &
    (df["reads_used"] <= q99_reads) &
    (df["covered_calls"] <= q99_calls) &
    (df["covered_bins"] <= q99_bins) &
    (df["promoter_covered_calls"] <= q99_prom)
)

filters = {
    "P7_QCpass_strict": strict,
    "P7_QCpass_very_strict": very_strict,
    "P7_QCpass_no01_reference": strict & (~df["cluster"].isin([0, 1])),
    "P7_QCpass_no016_reference": strict & (~df["cluster"].isin([0, 1, 6])),
}

rows = []
for name, mask in filters.items():
    sub = df[mask].copy()
    sub.to_csv(OUT / f"snTAPS-FL-1.{name}.global_per_cb.tsv", sep="\t", index=False)
    sub[["CB"]].to_csv(OUT / f"snTAPS-FL-1.{name}.cells.txt", sep="\t", index=False, header=False)

    row = {
        "filter": name,
        "n_cells": len(sub),
        "cluster0": int((sub["cluster"] == 0).sum()),
        "cluster1": int((sub["cluster"] == 1).sum()),
        "cluster6": int((sub["cluster"] == 6).sum()),
        "cluster2to7": int((~sub["cluster"].isin([0,1])).sum()),
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
summary.to_csv(OUT / "P7_QCpass_filter_summary.tsv", sep="\t", index=False)
print(summary.to_string(index=False))
