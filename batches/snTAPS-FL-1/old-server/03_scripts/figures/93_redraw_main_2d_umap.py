#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ========= input =========
UMAP_TSV = Path("08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv")
OUTDIR = Path("09_results/pig_frontal_main_figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

# old cluster -> new cluster name
RENAME_MAP = {
    2: "Cluster1",
    3: "Cluster2",
    4: "Cluster3",
    5: "Cluster4",
    6: "Cluster5",
    7: "Cluster6",
}

COLOR_MAP = {
    2: "#1f77b4",  # blue
    3: "#2ca02c",  # green
    4: "#9467bd",  # purple
    5: "#e377c2",  # pink
    6: "#bcbd22",  # olive
    7: "#17becf",  # cyan
}

df = pd.read_csv(UMAP_TSV, sep="\t")
df["cluster"] = pd.to_numeric(df["cluster"], errors="coerce").astype(int)

# keep cluster2-7 only
main_df = df[df["cluster"].isin(RENAME_MAP.keys())].copy()
main_df["main_cluster"] = main_df["cluster"].map(RENAME_MAP)

# save mapping / table
map_df = pd.DataFrame({
    "old_cluster": list(RENAME_MAP.keys()),
    "new_cluster": [RENAME_MAP[x] for x in RENAME_MAP]
})
map_df.to_csv(OUTDIR / "pig_frontal_cluster_rename_map.tsv", sep="\t", index=False)

main_df.to_csv(OUTDIR / "pig_frontal_main_umap_metadata.tsv", sep="\t", index=False)

count_df = (
    main_df.groupby(["cluster", "main_cluster"])
    .size()
    .reset_index(name="n_cells")
    .sort_values("cluster")
)
count_df.to_csv(OUTDIR / "pig_frontal_main_cluster_counts.tsv", sep="\t", index=False)

# plot
plt.figure(figsize=(8.0, 6.6))
for old_cl in sorted(RENAME_MAP.keys()):
    sub = main_df[main_df["cluster"] == old_cl]
    plt.scatter(
        sub["UMAP1"], sub["UMAP2"],
        s=8,
        c=COLOR_MAP[old_cl],
        linewidths=0,
        alpha=0.95
    )

# label at median position
for old_cl in sorted(RENAME_MAP.keys()):
    sub = main_df[main_df["cluster"] == old_cl]
    x = sub["UMAP1"].median()
    y = sub["UMAP2"].median()
    plt.text(
        x, y, RENAME_MAP[old_cl],
        fontsize=16, fontweight="bold",
        ha="center", va="center", color="black"
    )

plt.xlabel("UMAP1", fontsize=14)
plt.ylabel("UMAP2", fontsize=14)
plt.title("Pig frontal lobe methylation UMAP", fontsize=18)
plt.tight_layout()
plt.savefig(OUTDIR / "Fig1_pig_frontal_main_2D_UMAP.png", dpi=300)
plt.savefig(OUTDIR / "Fig1_pig_frontal_main_2D_UMAP.pdf")
plt.close()

print("[INFO] done")
print("[INFO] output:", OUTDIR)
