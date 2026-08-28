#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path

UMAP = Path("08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv")

OUTDIR = Path("08_qc_summary/final_single_sample/cell_filter_reads100k_q95")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(UMAP, sep="\t")

cell_col = None
for c in ["CB", "cell", "barcode", "cell_barcode"]:
    if c in df.columns:
        cell_col = c
        break
if cell_col is None:
    raise ValueError(f"Cannot find cell barcode column. columns={list(df.columns)}")

if "reads_used" not in df.columns:
    raise ValueError("reads_used not found in UMAP/QC table")

q = df["reads_used"].quantile([0.10, 0.50, 0.90, 0.95, 0.99])
q95 = q.loc[0.95]

mask = (df["reads_used"] >= 100000) & (df["reads_used"] <= q95)
keep = df.loc[mask].copy()
removed_low = df.loc[df["reads_used"] < 100000].copy()
removed_high = df.loc[df["reads_used"] > q95].copy()

keep[[cell_col]].to_csv(OUTDIR / "keep_reads100k_q95.cells.txt", sep="\t", index=False, header=False)
keep.to_csv(OUTDIR / "keep_reads100k_q95.umap_qc.tsv", sep="\t", index=False)
removed_low.to_csv(OUTDIR / "removed_reads_lt_100k.tsv", sep="\t", index=False)
removed_high.to_csv(OUTDIR / "removed_reads_gt_q95.tsv", sep="\t", index=False)

summary = pd.DataFrame([
    ["original", len(df), 1.0],
    ["removed_reads_lt_100k", len(removed_low), len(removed_low) / len(df)],
    ["removed_reads_gt_q95", len(removed_high), len(removed_high) / len(df)],
    ["keep_reads100k_q95", len(keep), len(keep) / len(df)],
], columns=["group", "n_cells", "fraction"])

summary.to_csv(OUTDIR / "cell_filter_summary.tsv", sep="\t", index=False)

cutoff = pd.DataFrame({
    "metric": ["reads_used"],
    "q10": [q.loc[0.10]],
    "median": [q.loc[0.50]],
    "q90": [q.loc[0.90]],
    "q95": [q.loc[0.95]],
    "q99": [q.loc[0.99]],
    "lower_cutoff_used": [100000],
    "upper_cutoff_used": [q95],
})
cutoff.to_csv(OUTDIR / "reads_filter_cutoffs.tsv", sep="\t", index=False)

print("[DONE]", OUTDIR)
print(summary)
print(cutoff)
