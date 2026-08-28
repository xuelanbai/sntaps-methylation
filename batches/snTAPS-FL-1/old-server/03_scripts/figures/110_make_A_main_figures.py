from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ========= input =========
UMAP_FILE = Path("08_qc_summary/final_single_sample/region_umap_paper/P7_no01_promoter_residual_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv")
GROUP_Z = Path("08_qc_summary/final_single_sample/marker_pig_frontal_lobe/A_no01_rerun_marker/A_no01_marker_group_activity_zscore.tsv")
PROMOTER_Z = Path("08_qc_summary/final_single_sample/marker_pig_frontal_lobe/A_no01_rerun_marker/A_no01_marker_promoter_activity_zscore.tsv")

OUTDIR = Path("09_results/figures_A_main")
OUTDIR.mkdir(parents=True, exist_ok=True)

# ========= cluster annotation =========
cluster_to_celltype = {
    0: "Inhibitory neurons",
    1: "Excitatory neurons I",
    2: "Vascular-associated cells",
    3: "Oligodendrocyte precursor cells",
    4: "Astrocytes",
    5: "Excitatory neurons II",
}

cluster_to_label = {
    0: "Inhibitory\nneurons",
    1: "Excitatory\nneurons I",
    2: "Vascular-associated\ncells",
    3: "Oligodendrocyte precursor\ncells",
    4: "Astrocytes",
    5: "Excitatory\nneurons II",
}

cluster_to_color = {
    0: "#377eb8",  # blue
    1: "#4daf4a",  # green
    2: "#35c0d3",  # cyan
    3: "#e78ac3",  # pink
    4: "#9c6bcc",  # purple
    5: "#bdbd22",  # olive
}

# 手动优化后的标签位置，避免挡图、超边框
label_pos = {
    0: (5.6, -0.4),   # Inhibitory
    1: (9.5, 8.8),    # Excitatory I
    2: (-7.2, 7.5),   # Vascular
    3: (0.9, 5.2),    # OPC
    4: (-5.0, 14.8),  # Astrocytes
    5: (2.9, 15.1),   # Excitatory II
}

def find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"None of these columns found: {candidates}; available={df.columns.tolist()}")

def plot_umap():
    print("[INFO] plot UMAP")
    df = pd.read_csv(UMAP_FILE, sep="\t")
    xcol = find_col(df, ["UMAP1", "umap_1", "umap1"])
    ycol = find_col(df, ["UMAP2", "umap_2", "umap2"])
    ccol = find_col(df, ["cluster", "Cluster"])

    df[ccol] = df[ccol].astype(int)
    df["cell_type"] = df[ccol].map(cluster_to_celltype)
    df["color"] = df[ccol].map(cluster_to_color)

    # 输出一份带细胞类型注释的表，后面备用
    out_tsv = OUTDIR / "Fig1_porcine_frontal_lobe_sc_methylation_UMAP_celltypes.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(9.2, 8.0), dpi=300)
    ax.scatter(df[xcol], df[ycol], c=df["color"], s=5, linewidths=0, alpha=0.9)

    for k in sorted(cluster_to_label):
        ax.text(
            label_pos[k][0], label_pos[k][1], cluster_to_label[k],
            fontsize=10.5, fontweight="bold",
            ha="center", va="center", color="black"
        )

    ax.set_title("Porcine frontal lobe single-cell methylation UMAP", fontsize=20, pad=16)
    ax.set_xlabel("UMAP1", fontsize=16)
    ax.set_ylabel("UMAP2", fontsize=16)
    ax.tick_params(labelsize=12)

    # 稍微留白，防止文字顶边
    xmin, xmax = df[xcol].min(), df[xcol].max()
    ymin, ymax = df[ycol].min(), df[ycol].max()
    ax.set_xlim(xmin - 1.0, xmax + 1.2)
    ax.set_ylim(ymin - 1.0, ymax + 1.0)

    plt.tight_layout()
    fig.savefig(OUTDIR / "Fig1_porcine_frontal_lobe_sc_methylation_UMAP_celltypes.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / "Fig1_porcine_frontal_lobe_sc_methylation_UMAP_celltypes.png", bbox_inches="tight")
    plt.close(fig)

def load_heatmap_df(path):
    df = pd.read_csv(path, sep="\t", index_col=0)

    rename_map = {
        "C0": "Inhibitory\nneurons",
        "C1": "Excitatory\nneurons I",
        "C2": "Vascular-associated\ncells",
        "C3": "Oligodendrocyte precursor\ncells",
        "C4": "Astrocytes",
        "C5": "Excitatory\nneurons II",
        0: "Inhibitory\nneurons",
        1: "Excitatory\nneurons I",
        2: "Vascular-associated\ncells",
        3: "Oligodendrocyte precursor\ncells",
        4: "Astrocytes",
        5: "Excitatory\nneurons II",
    }

    new_cols = []
    for c in df.columns:
        c2 = rename_map.get(c, rename_map.get(str(c), c))
        new_cols.append(c2)
    df.columns = new_cols
    return df

def plot_heatmap(df, title, outfile_prefix, figsize=(8.6, 5.8), y_fontsize=10):
    data = df.values.astype(float)

    fig, ax = plt.subplots(figsize=figsize, dpi=300)
    im = ax.imshow(data, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)

    ax.set_xticks(np.arange(df.shape[1]))
    ax.set_xticklabels(df.columns, fontsize=10)
    ax.set_yticks(np.arange(df.shape[0]))
    ax.set_yticklabels(df.index, fontsize=y_fontsize)

    ax.set_title(title, fontsize=18, pad=14)
    ax.set_xlabel("Cell type", fontsize=14)
    ax.set_ylabel("", fontsize=12)

    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Relative promoter hypomethylation activity\nrow Z-score", fontsize=12)
    cbar.ax.tick_params(labelsize=10)

    plt.tight_layout()
    fig.savefig(OUTDIR / f"{outfile_prefix}.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / f"{outfile_prefix}.png", bbox_inches="tight")
    plt.close(fig)

def main():
    if not UMAP_FILE.exists():
        raise FileNotFoundError(f"UMAP file not found: {UMAP_FILE}")
    plot_umap()

    if GROUP_Z.exists():
        print("[INFO] plot group heatmap")
        df_group = load_heatmap_df(GROUP_Z)
        plot_heatmap(
            df_group,
            "Porcine frontal lobe single-cell methylation: marker group activity",
            "Fig2_marker_group_activity",
            figsize=(8.4, 5.6),
            y_fontsize=10
        )
    else:
        print(f"[WARN] group heatmap input missing: {GROUP_Z}")

    if PROMOTER_Z.exists():
        print("[INFO] plot promoter heatmap")
        df_promoter = load_heatmap_df(PROMOTER_Z)
        plot_heatmap(
            df_promoter,
            "Porcine frontal lobe single-cell methylation: promoter hypomethylation marker activity",
            "Fig3_promoter_hypomethylation_marker_activity",
            figsize=(9.6, 15.0),
            y_fontsize=8
        )
    else:
        print(f"[WARN] promoter heatmap input missing: {PROMOTER_Z}")

    print("[INFO] done")
    print("[INFO] outputs in:", OUTDIR)

if __name__ == "__main__":
    main()
