#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

UMAP_TSV = Path("08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv")
OUTDIR = Path("09_results/pig_frontal_main_figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

# old cluster -> new cluster label
NEW_CLUSTER_MAP = {
    2: "Cluster1",
    3: "Cluster2",
    4: "Cluster3",
    5: "Cluster4",
    6: "Cluster5",
    7: "Cluster6",
}

# old cluster -> cell type
TYPE_MAP = {
    2: "Vascular",
    3: "Neuronal",
    4: "Microglia",
    5: "Oligodendrocyte",
    6: "OPC",
    7: "Inhibitory",
}

# old cluster -> plot label
LABEL_MAP = {
    k: f"{NEW_CLUSTER_MAP[k]}\n{TYPE_MAP[k]}"
    for k in NEW_CLUSTER_MAP
}

COLOR_MAP = {
    2: "#1f77b4",
    3: "#2ca02c",
    4: "#9467bd",
    5: "#e377c2",
    6: "#bcbd22",
    7: "#17becf",
}

# label position fine-tune
TEXT_OFFSET = {
    2: (-0.20,  0.15),
    3: ( 0.00,  0.30),
    4: ( 0.00, -0.10),
    5: ( 0.10,  0.18),
    6: (-0.10, -0.08),
    7: ( 0.00,  0.18),
}

df = pd.read_csv(UMAP_TSV, sep="\t")
df["cluster"] = pd.to_numeric(df["cluster"], errors="coerce").astype(int)

main = df[df["cluster"].isin(NEW_CLUSTER_MAP.keys())].copy()
main["new_cluster"] = main["cluster"].map(NEW_CLUSTER_MAP)
main["cell_type"] = main["cluster"].map(TYPE_MAP)

# 输出映射表和计数表
map_df = pd.DataFrame({
    "old_cluster": list(NEW_CLUSTER_MAP.keys()),
    "new_cluster": [NEW_CLUSTER_MAP[k] for k in NEW_CLUSTER_MAP],
    "cell_type": [TYPE_MAP[k] for k in NEW_CLUSTER_MAP],
})
map_df.to_csv(OUTDIR / "pig_frontal_cluster_rename_map.tsv", sep="\t", index=False)

count_df = (
    main.groupby(["cluster", "new_cluster", "cell_type"])
    .size()
    .reset_index(name="n_cells")
    .sort_values("cluster")
)
count_df.to_csv(OUTDIR / "pig_frontal_cluster_counts_main.tsv", sep="\t", index=False)

# plot
plt.figure(figsize=(8.6, 6.8))

for cl in sorted(NEW_CLUSTER_MAP.keys()):
    sub = main[main["cluster"] == cl]
    plt.scatter(
        sub["UMAP1"], sub["UMAP2"],
        s=8,
        c=COLOR_MAP[cl],
        linewidths=0,
        alpha=0.95
    )

for cl in sorted(NEW_CLUSTER_MAP.keys()):
    sub = main[main["cluster"] == cl]
    x = sub["UMAP1"].median() + TEXT_OFFSET[cl][0]
    y = sub["UMAP2"].median() + TEXT_OFFSET[cl][1]

    plt.text(
        x, y,
        LABEL_MAP[cl],
        fontsize=13,
        fontweight="bold",
        ha="center",
        va="center",
        color="black",
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.60,
            boxstyle="round,pad=0.18"
        )
    )

plt.xlabel("UMAP1", fontsize=15)
plt.ylabel("UMAP2", fontsize=15)
plt.title("Pig frontal lobe methylation UMAP", fontsize=20)
plt.tight_layout()

plt.savefig(OUTDIR / "Fig1_pig_frontal_methylation_UMAP_cleanlabel.png", dpi=300)
plt.savefig(OUTDIR / "Fig1_pig_frontal_methylation_UMAP_cleanlabel.pdf")
plt.close()

print("[INFO] wrote:")
print(OUTDIR / "Fig1_pig_frontal_methylation_UMAP_cleanlabel.png")
print(OUTDIR / "Fig1_pig_frontal_methylation_UMAP_cleanlabel.pdf")
print(OUTDIR / "pig_frontal_cluster_rename_map.tsv")
print(OUTDIR / "pig_frontal_cluster_counts_main.tsv")
