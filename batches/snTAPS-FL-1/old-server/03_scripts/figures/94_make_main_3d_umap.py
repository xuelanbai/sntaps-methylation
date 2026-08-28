#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import umap

# ========= set paths =========
META_TSV = Path("08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv")

# 把这里改成你找到的PCA输入文件
PCA_INPUT = Path("YOUR_PCA_INPUT.tsv")

OUTDIR = Path("09_results/pig_frontal_main_figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

RENAME_MAP = {
    2: "Cluster1",
    3: "Cluster2",
    4: "Cluster3",
    5: "Cluster4",
    6: "Cluster5",
    7: "Cluster6",
}

COLOR_MAP = {
    2: "#1f77b4",
    3: "#2ca02c",
    4: "#9467bd",
    5: "#e377c2",
    6: "#bcbd22",
    7: "#17becf",
}

# ========= read meta =========
meta = pd.read_csv(META_TSV, sep="\t")
meta["cluster"] = pd.to_numeric(meta["cluster"], errors="coerce").astype(int)
meta = meta[meta["cluster"].isin(RENAME_MAP.keys())].copy()
meta["main_cluster"] = meta["cluster"].map(RENAME_MAP)
meta["CB"] = meta["CB"].astype(str)

if not PCA_INPUT.exists():
    raise FileNotFoundError(f"PCA input not found: {PCA_INPUT}")

pca_df = pd.read_csv(PCA_INPUT, sep="\t")
if "CB" not in pca_df.columns:
    raise ValueError(f"PCA file must contain column 'CB'. Current columns: {pca_df.columns.tolist()}")

pca_df["CB"] = pca_df["CB"].astype(str)

# intersect cells
use = meta.merge(pca_df, on="CB", how="inner")
if use.shape[0] == 0:
    raise RuntimeError("No overlapping cells between metadata and PCA input.")

# select numeric PCA columns
feature_cols = []
for c in use.columns:
    if c in meta.columns:
        continue
    if c == "CB":
        continue
    if pd.api.types.is_numeric_dtype(use[c]):
        feature_cols.append(c)

if len(feature_cols) < 3:
    raise RuntimeError(f"Too few numeric columns for 3D UMAP: {feature_cols}")

X = use[feature_cols].to_numpy()

# ========= run 3D UMAP =========
reducer = umap.UMAP(
    n_components=3,
    n_neighbors=30,
    min_dist=0.3,
    metric="euclidean",
    random_state=42,
)
emb = reducer.fit_transform(X)

out = use[["CB", "cluster", "main_cluster"]].copy()
out["UMAP1"] = emb[:, 0]
out["UMAP2"] = emb[:, 1]
out["UMAP3"] = emb[:, 2]
out.to_csv(OUTDIR / "pig_frontal_main_3d_umap.tsv", sep="\t", index=False)

# ========= 3 views =========
def save_view(elev, azim, name):
    fig = plt.figure(figsize=(8.5, 7.2))
    ax = fig.add_subplot(111, projection='3d')

    for old_cl in sorted(RENAME_MAP.keys()):
        sub = out[out["cluster"] == old_cl]
        ax.scatter(
            sub["UMAP1"], sub["UMAP2"], sub["UMAP3"],
            s=8, alpha=0.95, c=COLOR_MAP[old_cl], depthshade=False
        )

    for old_cl in sorted(RENAME_MAP.keys()):
        sub = out[out["cluster"] == old_cl]
        x = sub["UMAP1"].median()
        y = sub["UMAP2"].median()
        z = sub["UMAP3"].median()
        ax.text(x, y, z, RENAME_MAP[old_cl], fontsize=12, fontweight='bold')

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_zlabel("UMAP3")
    ax.set_title("Pig frontal lobe methylation 3D UMAP", pad=20)
    ax.view_init(elev=elev, azim=azim)
    plt.tight_layout()
    plt.savefig(OUTDIR / name, dpi=300)
    plt.close()

save_view(20, -60, "Fig2_pig_frontal_main_3D_UMAP_view1.png")
save_view(15,  20, "Fig2_pig_frontal_main_3D_UMAP_view2.png")
save_view(35, 120, "Fig2_pig_frontal_main_3D_UMAP_view3.png")

# ========= optional interactive HTML =========
try:
    import plotly.express as px
    fig = px.scatter_3d(
        out,
        x="UMAP1", y="UMAP2", z="UMAP3",
        color="main_cluster",
        hover_data=["CB", "cluster"],
        title="Pig frontal lobe methylation 3D UMAP"
    )
    fig.write_html(str(OUTDIR / "Fig2_pig_frontal_main_3D_UMAP_interactive.html"))
    print("[INFO] interactive html saved")
except Exception as e:
    print("[WARN] plotly html not generated:", e)

print("[INFO] done")
print("[INFO] output:", OUTDIR)
