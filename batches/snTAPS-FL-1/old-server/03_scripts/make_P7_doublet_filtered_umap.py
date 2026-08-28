#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

UMAP = "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"
SUSPECT = "08_qc_summary/final_single_sample/P7_cluster01_doublet_check/suspected_doublet_like_cluster01_cells.tsv"
OUTDIR = Path("08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_doubletFiltered")
OUTDIR.mkdir(parents=True, exist_ok=True)

PREFIX = "snTAPS-FL-1.P7_promoter_residual_doubletFiltered"

def find_col(df, names):
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n in df.columns:
            return n
        if n.lower() in lower:
            return lower[n.lower()]
    raise ValueError(f"Cannot find any of {names}; columns={list(df.columns)}")

umap = pd.read_csv(UMAP, sep="\t")
sus = pd.read_csv(SUSPECT, sep="\t")

cb_col = find_col(umap, ["cell_barcode", "CB", "barcode", "cell"])
cl_col = find_col(umap, ["cluster", "orig_cluster", "leiden", "kmeans_cluster"])
x_col = find_col(umap, ["UMAP1", "umap1"])
y_col = find_col(umap, ["UMAP2", "umap2"])

sus_cb_col = find_col(sus, ["cell_barcode", "CB", "barcode", "cell"])

umap[cb_col] = umap[cb_col].astype(str)
sus[sus_cb_col] = sus[sus_cb_col].astype(str)

remove = set(sus[sus_cb_col])
before_n = len(umap)

filtered = umap[~umap[cb_col].isin(remove)].copy()
after_n = len(filtered)

filtered.to_csv(OUTDIR / f"{PREFIX}.umap.tsv", sep="\t", index=False)

summary = pd.DataFrame({
    "metric": ["before_n_cells", "removed_suspected_multiplet_like", "after_n_cells", "removed_fraction"],
    "value": [before_n, before_n - after_n, after_n, (before_n - after_n) / before_n]
})
summary.to_csv(OUTDIR / f"{PREFIX}.filter_summary.tsv", sep="\t", index=False)

counts = (
    filtered[cl_col]
    .value_counts()
    .rename_axis("cluster")
    .reset_index(name="n_cells")
    .sort_values("cluster")
)
counts["fraction"] = counts["n_cells"] / counts["n_cells"].sum()
counts.to_csv(OUTDIR / f"{PREFIX}.cluster_counts.tsv", sep="\t", index=False)

# cluster UMAP
clusters = sorted(pd.to_numeric(filtered[cl_col], errors="coerce").dropna().astype(int).unique())
base_colors = plt.cm.tab10(np.linspace(0, 1, max(10, len(clusters))))
cmap = ListedColormap(base_colors[:len(clusters)])
norm = BoundaryNorm(np.arange(min(clusters)-0.5, max(clusters)+1.5, 1), cmap.N)

fig, ax = plt.subplots(figsize=(7, 6))
sc = ax.scatter(
    filtered[x_col], filtered[y_col],
    c=pd.to_numeric(filtered[cl_col], errors="coerce"),
    s=1.2, cmap=cmap, norm=norm, linewidths=0, alpha=0.9
)
ax.set_title("Promoter residual UMAP")
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
cbar = plt.colorbar(sc, ax=ax, ticks=clusters)
cbar.set_label("Cluster")
fig.tight_layout()
fig.savefig(OUTDIR / f"{PREFIX}.umap_by_cluster.png", dpi=300)
fig.savefig(OUTDIR / f"{PREFIX}.umap_by_cluster.pdf")
plt.close(fig)

# density hexbin
fig, ax = plt.subplots(figsize=(7, 6))
hb = ax.hexbin(filtered[x_col], filtered[y_col], gridsize=70, mincnt=1)
ax.set_title("Promoter residual UMAP density")
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
cbar = plt.colorbar(hb, ax=ax)
cbar.set_label("Cell density")
fig.tight_layout()
fig.savefig(OUTDIR / f"{PREFIX}.umap_density_hexbin.png", dpi=300)
fig.savefig(OUTDIR / f"{PREFIX}.umap_density_hexbin.pdf")
plt.close(fig)

# optional QC colored panels if columns exist
qc_cols = [
    ("covered_bins", "Colored by Covered bins"),
    ("covered_calls", "Colored by Covered calls"),
    ("modified_rate", "Colored by Modified rate"),
    ("reads_used", "Colored by Reads used"),
]

for col, title in qc_cols:
    if col not in filtered.columns:
        continue
    vals = pd.to_numeric(filtered[col], errors="coerce")
    lo, hi = vals.quantile(0.01), vals.quantile(0.99)
    fig, ax = plt.subplots(figsize=(7, 6))
    sc = ax.scatter(
        filtered[x_col], filtered[y_col],
        c=vals.clip(lo, hi),
        s=1.2, cmap="viridis", linewidths=0, alpha=0.9
    )
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label(col)
    fig.tight_layout()
    fig.savefig(OUTDIR / f"{PREFIX}.umap_by_{col}.png", dpi=300)
    fig.savefig(OUTDIR / f"{PREFIX}.umap_by_{col}.pdf")
    plt.close(fig)

print("Done")
print(summary)
print("Output:", OUTDIR)
