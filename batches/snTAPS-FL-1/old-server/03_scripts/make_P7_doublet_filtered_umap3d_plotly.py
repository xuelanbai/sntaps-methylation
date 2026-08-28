#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

UMAP3D = Path("08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100_3D/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100_3D.umap.tsv")
KEEP2D = Path("08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_doubletFiltered/snTAPS-FL-1.P7_promoter_residual_doubletFiltered.umap.tsv")

OUTDIR = Path("08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100_3D_doubletFiltered")
OUTDIR.mkdir(parents=True, exist_ok=True)

PREFIX = "snTAPS-FL-1.P7_promoter_residual_3D_doubletFiltered"

def find_col(df, candidates):
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c in df.columns:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    raise ValueError(f"Cannot find any of {candidates}; columns={list(df.columns)}")

print("[INFO] reading 3D UMAP")
umap3d = pd.read_csv(UMAP3D, sep="\t")

print("[INFO] reading filtered 2D UMAP as keep-cell list")
keep2d = pd.read_csv(KEEP2D, sep="\t")

cb3 = find_col(umap3d, ["cell_barcode", "CB", "barcode", "cell", "cell_id"])
cb2 = find_col(keep2d, ["cell_barcode", "CB", "barcode", "cell", "cell_id"])
cl2 = find_col(keep2d, ["cluster", "orig_cluster", "leiden", "kmeans_cluster"])

umap3d[cb3] = umap3d[cb3].astype(str)
keep2d[cb2] = keep2d[cb2].astype(str)

# 用doubletFiltered结果中的细胞作为保留名单
keep_cells = set(keep2d[cb2])

before_n = len(umap3d)
filtered3d = umap3d[umap3d[cb3].isin(keep_cells)].copy()
after_n = len(filtered3d)

# 用2D doubletFiltered中的cluster标签覆盖3D中的cluster标签，保证和论文主图一致
cluster_map = keep2d.set_index(cb2)[cl2].to_dict()
filtered3d["cluster"] = filtered3d[cb3].map(cluster_map)

summary = pd.DataFrame({
    "metric": ["before_n_cells_3D", "after_n_cells_3D_doubletFiltered", "removed_n", "removed_fraction"],
    "value": [before_n, after_n, before_n - after_n, (before_n - after_n) / before_n]
})
summary.to_csv(OUTDIR / f"{PREFIX}.filter_summary.tsv", sep="\t", index=False)

filtered3d.to_csv(OUTDIR / f"{PREFIX}.umap.tsv", sep="\t", index=False)

counts = (
    filtered3d["cluster"]
    .value_counts()
    .rename_axis("cluster")
    .reset_index(name="n_cells")
    .sort_values("cluster")
)
counts["fraction"] = counts["n_cells"] / counts["n_cells"].sum()
counts.to_csv(OUTDIR / f"{PREFIX}.cluster_counts.tsv", sep="\t", index=False)

print("[INFO] before_n =", before_n)
print("[INFO] after_n =", after_n)
print("[INFO] removed =", before_n - after_n)
print("[INFO] output =", OUTDIR)

# Plotly HTML
try:
    import plotly.express as px
except ImportError:
    raise SystemExit(
        "\n[ERROR] plotly not installed in current environment.\n"
        "Try: python3 -m pip install --user plotly\n"
    )

plot_dir = OUTDIR / "plotly_html"
plot_dir.mkdir(exist_ok=True)

# cluster 3D
df = filtered3d.copy()
df["cluster"] = df["cluster"].astype(str)

fig = px.scatter_3d(
    df,
    x="UMAP1",
    y="UMAP2",
    z="UMAP3",
    color="cluster",
    title="Promoter residual 3D UMAP after suspected multiplet-like cell filtering",
    opacity=0.75,
    hover_data=[cb3, "cluster"] if cb3 in df.columns else ["cluster"]
)
fig.update_traces(marker=dict(size=2))
fig.update_layout(
    scene=dict(
        xaxis_title="UMAP1",
        yaxis_title="UMAP2",
        zaxis_title="UMAP3"
    ),
    legend_title_text="Cluster"
)
fig.write_html(plot_dir / f"{PREFIX}.by_cluster.html")

# QC colored 3D
for col in ["reads_used", "covered_bins", "covered_calls", "modified_rate"]:
    if col not in df.columns:
        continue
    vals = pd.to_numeric(df[col], errors="coerce")
    lo, hi = vals.quantile(0.01), vals.quantile(0.99)
    df[f"{col}_clip"] = vals.clip(lo, hi)

    fig = px.scatter_3d(
        df,
        x="UMAP1",
        y="UMAP2",
        z="UMAP3",
        color=f"{col}_clip",
        title=f"Promoter residual 3D UMAP colored by {col}",
        opacity=0.75,
        hover_data=[cb3, "cluster", col] if cb3 in df.columns else ["cluster", col],
        color_continuous_scale="Viridis"
    )
    fig.update_traces(marker=dict(size=2))
    fig.update_layout(
        scene=dict(
            xaxis_title="UMAP1",
            yaxis_title="UMAP2",
            zaxis_title="UMAP3"
        )
    )
    fig.write_html(plot_dir / f"{PREFIX}.by_{col}.html")

print("[INFO] plotly html output =", plot_dir)
