#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(".")
OUT = BASE / "09_results/P7_cluster2to7_initial_figures"
OUT.mkdir(parents=True, exist_ok=True)

P7_UMAP = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"

GLOBAL_QC = BASE / "08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.global_per_cb.tsv"

QCPASS_STRICT = BASE / "08_qc_summary/final_single_sample/cell_filter_P7_QCpass/snTAPS-FL-1.P7_QCpass_strict.global_per_cb.tsv"
QCPASS_NO01 = BASE / "08_qc_summary/final_single_sample/cell_filter_P7_QCpass/snTAPS-FL-1.P7_QCpass_no01_reference.global_per_cb.tsv"

MARKER_Z = BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/P7_marker_module_score/all_cluster_marker_trend/P7_QCpass_strict.all_clusters.weighted_activity.row_zscore.tsv"
MARKER_DELTA = BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/P7_marker_module_score/all_cluster_marker_trend/P7_QCpass_strict.all_clusters.module_delta_summary.tsv"
MARKER_ANNO = BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/P7_marker_module_score/all_cluster_marker_trend/P7_QCpass_strict.all_clusters.trend_annotation.tsv"

trend_label = {
    2: "C2 vascular/glial-like weak",
    3: "C3 neuronal-like/unclear",
    4: "C4 microglia-like",
    5: "C5 oligodendrocyte-like",
    6: "C6 OPC-like weak",
    7: "C7 inhibitory-like weak",
}

print("[INFO] read P7 UMAP", flush=True)
p7 = pd.read_csv(P7_UMAP, sep="\t")
p7["CB"] = p7["CB"].astype(str)
p7["cluster"] = pd.to_numeric(p7["cluster"], errors="coerce").astype(int)

p7_no01 = p7[~p7["cluster"].isin([0, 1])].copy()
p7_no01["trend_label"] = p7_no01["cluster"].map(trend_label)

p7.to_csv(OUT / "P7_all_cells.umap.tsv", sep="\t", index=False)
p7_no01.to_csv(OUT / "P7_cluster2to7_cells.umap.tsv", sep="\t", index=False)
p7_no01[["CB"]].to_csv(OUT / "P7_cluster2to7.CB_list.txt", sep="\t", index=False, header=False)

def plot_umap(df, outfile, title, annotate=False):
    plt.figure(figsize=(7, 6))
    sc = plt.scatter(
        df["UMAP1"],
        df["UMAP2"],
        c=df["cluster"],
        s=3,
        cmap="tab10",
        linewidths=0
    )
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(title)
    cbar = plt.colorbar(sc)
    cbar.set_label("P7 cluster")

    if annotate:
        for cl, sub in df.groupby("cluster"):
            x = sub["UMAP1"].median()
            y = sub["UMAP2"].median()
            plt.text(x, y, f"C{cl}", fontsize=10, weight="bold")

    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

plot_umap(
    p7,
    OUT / "Fig1A_P7_all_clusters_UMAP.png",
    "P7 original UMAP: all clusters",
    annotate=True
)

plot_umap(
    p7_no01,
    OUT / "Fig1B_P7_cluster2to7_UMAP_remove_cluster0_1.png",
    "P7 UMAP after removing cluster0/1",
    annotate=True
)

# cluster counts
counts_all = p7["cluster"].value_counts().sort_index().reset_index()
counts_all.columns = ["cluster", "n_cells"]
counts_all["fraction"] = counts_all["n_cells"] / counts_all["n_cells"].sum()
counts_all.to_csv(OUT / "P7_all_cluster_counts.tsv", sep="\t", index=False)

counts_no01 = p7_no01["cluster"].value_counts().sort_index().reset_index()
counts_no01.columns = ["cluster", "n_cells"]
counts_no01["fraction"] = counts_no01["n_cells"] / counts_no01["n_cells"].sum()
counts_no01.to_csv(OUT / "P7_cluster2to7_cluster_counts.tsv", sep="\t", index=False)

plt.figure(figsize=(6, 4))
plt.bar([f"C{x}" for x in counts_no01["cluster"]], counts_no01["n_cells"])
plt.xlabel("P7 cluster")
plt.ylabel("n cells")
plt.title("Cell numbers after removing cluster0/1")
for i, v in enumerate(counts_no01["n_cells"]):
    plt.text(i, v, str(v), ha="center", va="bottom", fontsize=8)
plt.tight_layout()
plt.savefig(OUT / "Fig1C_P7_cluster2to7_cell_counts.png", dpi=300)
plt.close()

# retention table and bar
rows = []
if GLOBAL_QC.exists():
    n_global = len(pd.read_csv(GLOBAL_QC, sep="\t", usecols=["CB"]))
else:
    n_global = np.nan

rows.append({"stage": "All CpG-covered CB", "n_cells": n_global})
rows.append({"stage": "strict2 QC", "n_cells": 73676})
rows.append({"stage": "P7 clustered", "n_cells": len(p7)})
rows.append({"stage": "P7 remove C0/C1", "n_cells": len(p7_no01)})

if QCPASS_STRICT.exists():
    rows.append({"stage": "P7 QCpass strict", "n_cells": len(pd.read_csv(QCPASS_STRICT, sep="\t", usecols=["CB"]))})
if QCPASS_NO01.exists():
    rows.append({"stage": "P7 QCpass no C0/C1", "n_cells": len(pd.read_csv(QCPASS_NO01, sep="\t", usecols=["CB"]))})

ret = pd.DataFrame(rows)
ret.to_csv(OUT / "P7_cell_retention_summary.tsv", sep="\t", index=False)

plt.figure(figsize=(8, 4))
plt.bar(range(len(ret)), ret["n_cells"])
plt.xticks(range(len(ret)), ret["stage"], rotation=35, ha="right")
plt.ylabel("n CB/cells")
plt.title("Cell/CB retention across filtering steps")
for i, v in enumerate(ret["n_cells"]):
    if pd.notna(v):
        plt.text(i, v, f"{int(v):,}", ha="center", va="bottom", fontsize=8)
plt.tight_layout()
plt.savefig(OUT / "Fig1D_cell_retention_bar.png", dpi=300)
plt.close()

# QC boxplots
qc_metrics = [
    "reads_used",
    "covered_calls",
    "covered_bins",
    "modified_rate",
]
qc_metrics = [x for x in qc_metrics if x in p7.columns]

def boxplot_by_cluster(df, metrics, prefix, title_prefix):
    clusters = sorted(df["cluster"].unique())
    for m in metrics:
        data = [pd.to_numeric(df.loc[df["cluster"] == cl, m], errors="coerce").dropna().values for cl in clusters]
        plt.figure(figsize=(7, 4.5))
        plt.boxplot(data, labels=[f"C{x}" for x in clusters], showfliers=False)
        plt.xlabel("P7 cluster")
        plt.ylabel(m)
        plt.title(f"{title_prefix}: {m}")
        plt.tight_layout()
        plt.savefig(OUT / f"{prefix}_{m}_boxplot.png", dpi=300)
        plt.close()

boxplot_by_cluster(
    p7,
    qc_metrics,
    "Fig2A_P7_all_clusters_QC",
    "P7 all clusters QC"
)

boxplot_by_cluster(
    p7_no01,
    qc_metrics,
    "Fig2B_P7_cluster2to7_QC",
    "P7 cluster2-7 QC"
)

# Marker module heatmap for cluster2-7
if MARKER_Z.exists():
    z = pd.read_csv(MARKER_Z, sep="\t", index_col=0)

    # columns may be 0,1,2 or C0,C1...
    rename = {}
    for c in z.columns:
        s = str(c)
        if s.startswith("C"):
            rename[c] = s
        else:
            rename[c] = f"C{s}"
    z = z.rename(columns=rename)

    keep_cols = [f"C{x}" for x in [2, 3, 4, 5, 6, 7] if f"C{x}" in z.columns]
    z2 = z[keep_cols].copy()

    if MARKER_DELTA.exists():
        delta = pd.read_csv(MARKER_DELTA, sep="\t")
        order = delta["group"].tolist()
        order = [x for x in order if x in z2.index]
        z2 = z2.loc[order]

    z2.to_csv(OUT / "P7_cluster2to7_marker_module_row_zscore.tsv", sep="\t")

    plt.figure(figsize=(7, max(4, 0.45 * z2.shape[0])))
    im = plt.imshow(z2.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    plt.xticks(range(z2.shape[1]), z2.columns, rotation=45, ha="right")
    plt.yticks(range(z2.shape[0]), z2.index)
    plt.xlabel("P7 cluster")
    plt.ylabel("Pig frontal lobe marker module")
    plt.title("Marker module trend in P7 cluster2-7")
    plt.colorbar(im, label="row z-score")
    plt.tight_layout()
    plt.savefig(OUT / "Fig2C_P7_cluster2to7_marker_module_heatmap.png", dpi=300)
    plt.close()

if MARKER_DELTA.exists():
    delta = pd.read_csv(MARKER_DELTA, sep="\t")
    delta.to_csv(OUT / "P7_marker_module_delta_summary.tsv", sep="\t", index=False)

    plt.figure(figsize=(7, 4.5))
    d = delta.sort_values("delta", ascending=True)
    plt.barh(d["group"], d["delta"])
    plt.xlabel("max-min weighted activity across clusters")
    plt.ylabel("marker module")
    plt.title("Marker module separation strength")
    plt.tight_layout()
    plt.savefig(OUT / "Fig2D_marker_module_delta_bar.png", dpi=300)
    plt.close()

if MARKER_ANNO.exists():
    anno = pd.read_csv(MARKER_ANNO, sep="\t")
    anno.to_csv(OUT / "P7_all_cluster_marker_trend_annotation.tsv", sep="\t", index=False)

# simplified annotation table for cluster2-7
anno2 = pd.DataFrame([
    {"cluster": 2, "trend_annotation": "vascular/glial-like weak", "confidence": "weak"},
    {"cluster": 3, "trend_annotation": "neuronal-like / unclear", "confidence": "very weak"},
    {"cluster": 4, "trend_annotation": "microglia-like", "confidence": "moderate exploratory"},
    {"cluster": 5, "trend_annotation": "oligodendrocyte-like", "confidence": "weak"},
    {"cluster": 6, "trend_annotation": "OPC-like / upper-layer-like weak", "confidence": "very weak"},
    {"cluster": 7, "trend_annotation": "inhibitory-neuron-like weak", "confidence": "very weak"},
])
anno2.to_csv(OUT / "P7_cluster2to7_trend_annotation_for_draft.tsv", sep="\t", index=False)

print("[INFO] wrote figures and tables to:", OUT)
print("[INFO] P7 total cells:", len(p7))
print("[INFO] P7 cluster2-7 cells:", len(p7_no01))
print("[INFO] outputs:")
for f in sorted(OUT.iterdir()):
    print("  ", f)
