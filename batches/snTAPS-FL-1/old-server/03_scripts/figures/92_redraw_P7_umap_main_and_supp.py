#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(".")
OUT = BASE / "09_results/P7_cluster2to7_initial_figures"
OUT.mkdir(parents=True, exist_ok=True)

UMAP = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"

# 主图：old2-7 -> new Cluster1-6
RENAME_MAP = {
    2: "Cluster1",
    3: "Cluster2",
    4: "Cluster3",
    5: "Cluster4",
    6: "Cluster5",
    7: "Cluster6",
}

# 颜色尽量简洁、区分度够，但不用右侧colorbar
COLOR_MAP = {
    2: "#1f77b4",  # blue
    3: "#2ca02c",  # green
    4: "#9467bd",  # purple
    5: "#e377c2",  # pink
    6: "#bcbd22",  # olive
    7: "#17becf",  # cyan
}

SUPP_COLOR_MAP = {
    0: "#bdbdbd",
    1: "#969696",
    2: "#1f77b4",
    3: "#2ca02c",
    4: "#9467bd",
    5: "#e377c2",
    6: "#bcbd22",
    7: "#17becf",
}

df = pd.read_csv(UMAP, sep="\t")
df["CB"] = df["CB"].astype(str)
df["cluster"] = pd.to_numeric(df["cluster"], errors="coerce").astype(int)

# 输出旧->新编号对照表
map_df = pd.DataFrame({
    "old_cluster": list(RENAME_MAP.keys()),
    "new_cluster": [RENAME_MAP[k] for k in RENAME_MAP.keys()]
})
map_df.to_csv(OUT / "P7_main_cluster_rename_map.tsv", sep="\t", index=False)

# ========= 主图：仅cluster2-7 =========
main_df = df[df["cluster"].isin(RENAME_MAP.keys())].copy()
main_df["new_cluster"] = main_df["cluster"].map(RENAME_MAP)
main_df.to_csv(OUT / "P7_cluster2to7_renamed_umap.tsv", sep="\t", index=False)

plt.figure(figsize=(7.2, 6.0))

for old_cl in sorted(RENAME_MAP.keys()):
    sub = main_df[main_df["cluster"] == old_cl]
    plt.scatter(
        sub["UMAP1"], sub["UMAP2"],
        s=7,
        c=COLOR_MAP[old_cl],
        linewidths=0,
        alpha=0.95
    )

# 直接标注，不放legend/colorbar
for old_cl in sorted(RENAME_MAP.keys()):
    sub = main_df[main_df["cluster"] == old_cl]
    x = sub["UMAP1"].median()
    y = sub["UMAP2"].median()
    plt.text(
        x, y,
        RENAME_MAP[old_cl],
        fontsize=12,
        fontweight="bold",
        ha="center",
        va="center",
        color="black"
    )

plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title("P7 UMAP after removing cluster0/1")
plt.tight_layout()
plt.savefig(OUT / "FigMain_P7_cluster2to7_renamed_umap.png", dpi=300)
plt.savefig(OUT / "FigMain_P7_cluster2to7_renamed_umap.pdf")
plt.close()

# ========= 附录图：原始0-7 =========
supp_df = df.copy()

plt.figure(figsize=(7.2, 6.0))

for cl in sorted(supp_df["cluster"].unique()):
    sub = supp_df[supp_df["cluster"] == cl]
    plt.scatter(
        sub["UMAP1"], sub["UMAP2"],
        s=5,
        c=SUPP_COLOR_MAP.get(cl, "#333333"),
        linewidths=0,
        alpha=0.85
    )

for cl in sorted(supp_df["cluster"].unique()):
    sub = supp_df[supp_df["cluster"] == cl]
    x = sub["UMAP1"].median()
    y = sub["UMAP2"].median()
    plt.text(
        x, y,
        f"C{cl}",
        fontsize=11,
        fontweight="bold",
        ha="center",
        va="center",
        color="black"
    )

plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title("P7 original UMAP with old cluster IDs")
plt.tight_layout()
plt.savefig(OUT / "FigSup_P7_allcluster0to7_oldID_umap.png", dpi=300)
plt.savefig(OUT / "FigSup_P7_allcluster0to7_oldID_umap.pdf")
plt.close()

# ========= 再做一个主图细胞数表 =========
count_df = main_df.groupby(["cluster", "new_cluster"]).size().reset_index(name="n_cells")
count_df = count_df.sort_values("cluster")
count_df.to_csv(OUT / "P7_cluster2to7_renamed_counts.tsv", sep="\t", index=False)

print("[INFO] done")
print("[INFO] main figure:", OUT / "FigMain_P7_cluster2to7_renamed_umap.png")
print("[INFO] supplementary figure:", OUT / "FigSup_P7_allcluster0to7_oldID_umap.png")
print("[INFO] rename map:", OUT / "P7_main_cluster_rename_map.tsv")
print("[INFO] renamed counts:", OUT / "P7_cluster2to7_renamed_counts.tsv")
