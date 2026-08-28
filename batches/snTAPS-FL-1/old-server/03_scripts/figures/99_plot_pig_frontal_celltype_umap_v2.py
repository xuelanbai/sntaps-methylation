#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ===== input / output =====
INFILE = Path(
    "08_qc_summary/final_single_sample/region_umap_paper/"
    "P7_no01_promoter_residual_dropPC1_smooth100/"
    "snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"
)

OUTDIR = Path("09_results/figures_main_pig_frontal")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_PDF = OUTDIR / "Fig1_pig_frontal_celltype_UMAP.pdf"
OUT_PNG = OUTDIR / "Fig1_pig_frontal_celltype_UMAP.png"

# ===== title =====
TITLE = "Porcine frontal lobe single-cell methylation UMAP"

# ===== cluster -> cell type mapping =====
# 这里按你们正文想展示的主线结果来改
# 如果当前输入文件里的cluster不是1~6，而是0~5，就把左边数字改掉
CLUSTER_TO_CELLTYPE = {
    0: "Vascular cells",
    1: "Neurons",
    2: "Microglia",
    3: "Oligodendrocytes",
    4: "Oligodendrocyte precursor cells",
    5: "Inhibitory neurons",
}

# ===== colors =====
CELLTYPE_COLORS = {
    "Vascular cells": "#2C7FB8",
    "Neurons": "#33A02C",
    "Microglia": "#9467BD",
    "Oligodendrocytes": "#E377C2",
    "Oligodendrocyte precursor cells": "#BCBD22",
    "Inhibitory neurons": "#17BECF",
}

# ===== optional manual label positions =====
# 没写的会自动用中位数位置
LABEL_POS = {
    "Vascular cells": (-3.9, 8.5),
    "Neurons": (1.1, 7.0),
    "Microglia": (0.2, -0.4),
    "Oligodendrocytes": (10.3, 8.9),
    "Oligodendrocyte precursor cells": (-0.5, 3.3),
    "Inhibitory neurons": (7.0, 8.7),
}

def main():
    if not INFILE.exists():
        raise FileNotFoundError(f"UMAP file not found: {INFILE}")

    df = pd.read_csv(INFILE, sep="\t")
    required = {"UMAP1", "UMAP2", "cluster"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in input: {sorted(missing)}")

    # cluster列统一转成int
    df["cluster"] = pd.to_numeric(df["cluster"], errors="coerce")
    df = df.dropna(subset=["cluster", "UMAP1", "UMAP2"]).copy()
    df["cluster"] = df["cluster"].astype(int)

    # 只保留要展示的cluster
    df = df[df["cluster"].isin(CLUSTER_TO_CELLTYPE.keys())].copy()

    if df.empty:
        raise RuntimeError("No cells left after filtering target clusters.")

    df["celltype"] = df["cluster"].map(CLUSTER_TO_CELLTYPE)

    # 作图
    fig, ax = plt.subplots(figsize=(11, 8.5))

    plot_order = list(CLUSTER_TO_CELLTYPE.values())

    for celltype in plot_order:
        sub = df[df["celltype"] == celltype]
        if sub.empty:
            continue
        ax.scatter(
            sub["UMAP1"],
            sub["UMAP2"],
            s=6,
            alpha=0.85,
            linewidths=0,
            label=celltype,
            color=CELLTYPE_COLORS.get(celltype, "gray"),
            rasterized=True,
        )

    # 标签
    for celltype in plot_order:
        sub = df[df["celltype"] == celltype]
        if sub.empty:
            continue

        if celltype in LABEL_POS:
            x0, y0 = LABEL_POS[celltype]
        else:
            x0 = float(np.median(sub["UMAP1"]))
            y0 = float(np.median(sub["UMAP2"]))

        ax.text(
            x0,
            y0,
            celltype,
            fontsize=16,
            fontweight="bold",
            ha="center",
            va="center",
            color="black",
            bbox=dict(
                boxstyle="round,pad=0.18",
                facecolor="white",
                edgecolor="none",
                alpha=0.45,
            ),
        )

    ax.set_title(TITLE, fontsize=22, pad=18)
    ax.set_xlabel("UMAP1", fontsize=18)
    ax.set_ylabel("UMAP2", fontsize=18)
    ax.tick_params(axis="both", labelsize=13)

    # 图例放右侧，避免挡住图
    handles, labels = ax.get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = h

    ax.legend(
        uniq.values(),
        uniq.keys(),
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=14,
        markerscale=2.5,
        handletextpad=0.5,
        borderaxespad=0.6,
    )

    plt.tight_layout()
    fig.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[INFO] saved pdf: {OUT_PDF}")
    print(f"[INFO] saved png: {OUT_PNG}")
    print("[INFO] cluster counts:")
    print(df["cluster"].value_counts().sort_index().to_string())
    print("[INFO] cell type counts:")
    print(df["celltype"].value_counts().to_string())

if __name__ == "__main__":
    main()