#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(".")
UMAP = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"
QC = BASE / "08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.global_per_cb.tsv"
OUT = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_realcell_filter_diagnosis"
OUT.mkdir(parents=True, exist_ok=True)

u = pd.read_csv(UMAP, sep="\t")
q = pd.read_csv(QC, sep="\t")
u["CB"] = u["CB"].astype(str)
q["CB"] = q["CB"].astype(str)

# 合并QC列
qc_cols = ["reads_used", "modified_calls", "covered_calls", "modified_rate", "covered_bins"]
missing = [c for c in qc_cols if c not in u.columns and c in q.columns]
if missing:
    u = u.merge(q[["CB"] + missing], on="CB", how="left")

for c in qc_cols:
    if c in u.columns:
        u[c] = pd.to_numeric(u[c], errors="coerce")

u["cluster"] = pd.to_numeric(u["cluster"], errors="coerce").astype("Int64")

def plot_subset(df, name):
    if len(df) == 0:
        return

    plt.figure(figsize=(7,6))
    sc = plt.scatter(df["UMAP1"], df["UMAP2"], c=df["cluster"], s=3, cmap="tab10")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(f"{name} colored by old P7 cluster")
    plt.colorbar(sc, label="cluster")
    plt.tight_layout()
    plt.savefig(OUT / f"{name}.oldP7_by_cluster.png", dpi=220)
    plt.close()

    for col in ["reads_used", "covered_calls", "covered_bins", "modified_rate"]:
        if col in df.columns:
            plt.figure(figsize=(7,6))
            sc = plt.scatter(df["UMAP1"], df["UMAP2"], c=df[col], s=3)
            plt.xlabel("UMAP1")
            plt.ylabel("UMAP2")
            plt.title(f"{name} colored by {col}")
            plt.colorbar(sc, label=col)
            plt.tight_layout()
            plt.savefig(OUT / f"{name}.oldP7_by_{col}.png", dpi=220)
            plt.close()

    df["cluster"].value_counts().sort_index().to_csv(
        OUT / f"{name}.cluster_counts.tsv",
        sep="\t",
        header=["n_cells"]
    )

    df.to_csv(OUT / f"{name}.oldP7_subset.tsv", sep="\t", index=False)

def summarize(df, name):
    row = {
        "filter": name,
        "n_cells": len(df),
        "cluster0": int((df["cluster"] == 0).sum()),
        "cluster1": int((df["cluster"] == 1).sum()),
    }
    row["cluster0_1"] = row["cluster0"] + row["cluster1"]
    row["cluster0_1_frac"] = row["cluster0_1"] / row["n_cells"] if row["n_cells"] else np.nan
    for c in ["reads_used", "covered_calls", "covered_bins", "modified_rate"]:
        if c in df.columns:
            row[f"{c}_median"] = df[c].median()
    return row

rows = []

# all P7
plot_subset(u, "P7_all")
rows.append(summarize(u, "P7_all"))

# reads >= 200k
f = u[u["reads_used"] >= 200000].copy()
plot_subset(f, "P7_reads_ge200k")
rows.append(summarize(f, "P7_reads_ge200k"))

# reads >= 250k / 300k也顺手看一下，不重跑UMAP
for cutoff in [250000, 300000]:
    f = u[u["reads_used"] >= cutoff].copy()
    plot_subset(f, f"P7_reads_ge{cutoff//1000}k")
    rows.append(summarize(f, f"P7_reads_ge{cutoff//1000}k"))

# expected_realcell文件
candidates = sorted(Path("08_qc_summary/final_single_sample").rglob("*expected_realcell_top*global_per_cb.tsv"))
print("[INFO] expected_realcell candidates:")
for p in candidates:
    print("  ", p)

for p in candidates:
    name = p.name.replace(".global_per_cb.tsv", "")
    cells = pd.read_csv(p, sep="\t", usecols=["CB"])
    cells["CB"] = cells["CB"].astype(str)
    df = u[u["CB"].isin(set(cells["CB"]))].copy()
    plot_subset(df, name)
    rows.append(summarize(df, name))

summary = pd.DataFrame(rows)
summary.to_csv(OUT / "P7_filter_summary.tsv", sep="\t", index=False)

print("[INFO] wrote", OUT)
print(summary)
