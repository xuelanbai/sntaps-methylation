#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN = Path("08_qc_summary/final_single_sample/cell_filter_reads100k_q95/keep_reads100k_q95.umap_qc.tsv")
OUTDIR = Path("08_qc_summary/final_single_sample/cell_filter_reads100k_q95/figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(IN, sep="\t")

def find_numeric_col(candidates):
    for c in candidates:
        if c in df.columns and pd.api.types.is_numeric_dtype(df[c]):
            return c
    return None

reads_col = find_numeric_col(["reads_used", "reads", "read_count", "n_reads"])
cytosine_col = find_numeric_col(["covered_calls", "cytosine_number", "cytosine_count", "C_number", "n_cytosine"])
bin_col = find_numeric_col(["covered_bins", "bin_number", "bins_number", "n_bins", "n_100kb_bins"])

# 自动寻找CpG number相关列
cpg_candidates = []
for c in df.columns:
    cl = c.lower()
    if "cpg" in cl and pd.api.types.is_numeric_dtype(df[c]):
        cpg_candidates.append(c)

# 避免优先选rate/ratio
cpg_number_col = None
for c in cpg_candidates:
    cl = c.lower()
    if not any(x in cl for x in ["rate", "ratio", "mean", "frac", "level"]):
        cpg_number_col = c
        break

metrics = []
if reads_col:
    metrics.append(("reads_per_cell", reads_col, "log10(reads per cell + 1)", True))
if cytosine_col:
    metrics.append(("cytosine_number_per_cell", cytosine_col, "log10(cytosine calls per cell + 1)", True))
if cpg_number_col:
    metrics.append(("cpg_number_per_cell", cpg_number_col, "log10(CpG number per cell + 1)", True))
if bin_col:
    metrics.append(("covered_100kb_bins_per_cell", bin_col, "100-kb bins per cell", False))

summary_rows = []

for name, col, xlabel, log_transform in metrics:
    x = df[col].dropna()
    x = x[x >= 0]
    if log_transform:
        y = np.log10(x + 1)
    else:
        y = x

    qs = y.quantile([0.10, 0.50, 0.90, 0.95])
    raw_qs = x.quantile([0.10, 0.50, 0.90, 0.95])

    summary_rows.append({
        "metric": name,
        "column": col,
        "n_cells": len(x),
        "raw_q10": raw_qs.loc[0.10],
        "raw_median": raw_qs.loc[0.50],
        "raw_q90": raw_qs.loc[0.90],
        "raw_q95": raw_qs.loc[0.95],
        "plot_q10": qs.loc[0.10],
        "plot_median": qs.loc[0.50],
        "plot_q90": qs.loc[0.90],
        "plot_q95": qs.loc[0.95],
    })

    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    ax.hist(y, bins=80, density=True, alpha=0.75)

    for qname, qv in [("q10", qs.loc[0.10]), ("median", qs.loc[0.50]), ("q90", qs.loc[0.90]), ("q95", qs.loc[0.95])]:
        ax.axvline(qv, linestyle="--", linewidth=1.2)
        ymax = ax.get_ylim()[1]
        ax.text(qv, ymax * 0.95, qname, rotation=90, va="top", ha="left", fontsize=8)

    ax.set_title(name.replace("_", " "))
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    fig.tight_layout()
    fig.savefig(OUTDIR / f"{name}.density.png")
    plt.close(fig)

summary = pd.DataFrame(summary_rows)
summary.to_csv(OUTDIR / "filtered_cell_qc_metric_summary.tsv", sep="\t", index=False)

with open(OUTDIR / "column_detection_report.txt", "w") as f:
    f.write(f"input={IN}\n")
    f.write(f"reads_col={reads_col}\n")
    f.write(f"cytosine_col={cytosine_col}\n")
    f.write(f"cpg_number_col={cpg_number_col}\n")
    f.write(f"bin_col={bin_col}\n")
    f.write("cpg_candidates=" + ",".join(cpg_candidates) + "\n")
    if cpg_number_col is None:
        f.write("WARNING: No CpG number column was detected in this table. Need to compute CpG number from CpG conversion/call table if required.\n")

print("[DONE]", OUTDIR)
print(summary)
print("[COLUMN REPORT]", OUTDIR / "column_detection_report.txt")
