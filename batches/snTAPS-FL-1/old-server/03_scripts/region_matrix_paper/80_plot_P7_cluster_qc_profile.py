#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(".")
IN = BASE / "08_qc_summary/final_single_sample/cell_filter_P7_promoter_coverage/snTAPS-FL-1.P7_with_promoter_coverage_qc.tsv"
OUT = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_cluster_qc_profile/plots"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(IN, sep="\t")
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

clusters = sorted(df["cluster"].dropna().unique())

for c in metrics:
    plot_df = [df.loc[df["cluster"] == cl, c].dropna().values for cl in clusters]

    plt.figure(figsize=(8, 5))
    plt.boxplot(plot_df, labels=[str(x) for x in clusters], showfliers=False)
    plt.xlabel("P7 cluster")
    plt.ylabel(c)
    plt.title(f"P7 cluster QC: {c}")
    plt.tight_layout()
    plt.savefig(OUT / f"P7_cluster_boxplot_{c}.png", dpi=220)
    plt.close()

print("[INFO] wrote plots to", OUT)
