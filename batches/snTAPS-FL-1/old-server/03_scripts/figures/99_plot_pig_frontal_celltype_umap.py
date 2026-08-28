import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ===== input =====
INFILE = "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"
OUTDIR = "09_results/figures_main_pig_frontal"
os.makedirs(OUTDIR, exist_ok=True)

# ===== read =====
df = pd.read_csv(INFILE, sep="\t")

# 需要包含：UMAP1, UMAP2, cluster_new
# 如果没有cluster_new，自己从cluster映射
if "cluster_new" not in df.columns:
    if "cluster" not in df.columns:
        raise ValueError(f"Need cluster_new or cluster column. Current columns: {df.columns.tolist()}")
    old2new = {2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6}
    df = df[df["cluster"].isin(old2new.keys())].copy()
    df["cluster_new"] = df["cluster"].map(old2new)

type_map = {
    1: "Vascular",
    2: "Neuronal",
    3: "Microglia",
    4: "Oligodendrocyte",
    5: "OPC",
    6: "Inhibitory"
}

df["cell_type"] = df["cluster_new"].map(type_map)

# ===== color palette =====
color_map = {
    "Vascular": "#1f77b4",
    "Neuronal": "#2ca02c",
    "Microglia": "#9467bd",
    "Oligodendrocyte": "#e377c2",
    "OPC": "#bcbd22",
    "Inhibitory": "#17becf",
}

order = ["Vascular", "Neuronal", "Microglia", "Oligodendrocyte", "OPC", "Inhibitory"]

# ===== label position =====
# 先用中位数作为默认中心，后面可手动微调
centers = (
    df.groupby("cell_type")[["UMAP1", "UMAP2"]]
    .median()
    .reset_index()
)

# 如果需要更美观，可以手动调位置
manual_pos = {
    "Vascular": (-3.7, 8.4),
    "Neuronal": (1.6, 7.1),
    "Microglia": (0.1, -0.4),
    "Oligodendrocyte": (11.3, 8.9),
    "OPC": (-0.5, 3.3),
    "Inhibitory": (6.9, 8.6),
}

# ===== plot =====
fig, ax = plt.subplots(figsize=(10.5, 7.2))

for ct in order:
    sub = df[df["cell_type"] == ct]
    ax.scatter(
        sub["UMAP1"],
        sub["UMAP2"],
        s=4,
        alpha=0.85,
        c=color_map[ct],
        label=ct,
        linewidths=0
    )

# 图中标签
for ct in order:
    if ct in manual_pos:
        x, y = manual_pos[ct]
    else:
        row = centers[centers["cell_type"] == ct].iloc[0]
        x, y = row["UMAP1"], row["UMAP2"]

    ax.text(
        x, y, ct,
        fontsize=14,
        fontweight="bold",
        ha="center",
        va="center",
        color="black"
    )

# 样式
ax.set_title("Pig frontal lobe methylation UMAP", fontsize=20, pad=12)
ax.set_xlabel("UMAP1", fontsize=14)
ax.set_ylabel("UMAP2", fontsize=14)
ax.tick_params(axis="both", labelsize=11)

# 去掉上右边框，更清爽
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# 图例放右边
leg = ax.legend(
    title=None,
    loc="center left",
    bbox_to_anchor=(1.02, 0.3),
    frameon=False,
    fontsize=14,
    markerscale=1.4,
    handletextpad=0.5,
    labelspacing=0.5,
    borderaxespad=0.0
)

plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "Fig1_pig_frontal_celltype_UMAP_pretty.pdf"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(OUTDIR, "Fig1_pig_frontal_celltype_UMAP_pretty.png"), dpi=300, bbox_inches="tight")
plt.close()

print("[INFO] done")