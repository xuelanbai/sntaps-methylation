import scanpy as sc
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base = Path("/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation")
indir = base / "08_qc_summary/final_single_sample/region_umap_leiden/promoter_cov10_leiden"
h5ad = indir / "snTAPS-FL-1.promoter_cov10.leiden.h5ad"

outdir = base / "08_qc_summary/final_single_sample/region_umap_leiden/promoter_cov10_leiden_coarse"
outdir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(h5ad)

qc_cols = ["reads_used", "covered_calls", "covered_bins", "modified_rate"]

settings = [
    (150, 0.02),
    (150, 0.05),
    (300, 0.02),
    (300, 0.05),
    (300, 0.10),
    (500, 0.02),
    (500, 0.05),
    (500, 0.10),
]

summary = []

for nn, res in settings:
    tag = f"nn{nn}_res{res}"

    sc.pp.neighbors(adata, n_neighbors=nn, use_rep="X", key_added=f"neighbors_{tag}")
    sc.tl.umap(
        adata,
        neighbors_key=f"neighbors_{tag}",
        min_dist=0.6,
        spread=1.2,
        random_state=1
    )
    sc.tl.leiden(
        adata,
        resolution=res,
        neighbors_key=f"neighbors_{tag}",
        key_added=f"leiden_{tag}"
    )

    emb = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs_names)
    emb.insert(0, "CB", emb.index)

    for c in qc_cols:
        emb[c] = adata.obs[c].values

    emb["cluster"] = adata.obs[f"leiden_{tag}"].astype(str).values
    emb.to_csv(outdir / f"snTAPS-FL-1.promoter_cov10.coarse_{tag}.umap.tsv", sep="\t", index=False)

    counts = emb["cluster"].value_counts().sort_index()
    counts.to_csv(outdir / f"snTAPS-FL-1.promoter_cov10.coarse_{tag}.cluster_counts.tsv", sep="\t")

    n_cluster = counts.shape[0]
    max_frac = counts.max() / counts.sum()
    small_frac = counts[counts < 500].sum() / counts.sum()
    summary.append([tag, nn, res, len(emb), n_cluster, counts.max(), max_frac, small_frac])

    corr = emb[["UMAP1", "UMAP2"] + qc_cols].corr(method="spearman")
    corr.to_csv(outdir / f"snTAPS-FL-1.promoter_cov10.coarse_{tag}.qc_spearman_corr.tsv", sep="\t")

    # cluster图
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    cats = sorted(emb["cluster"].unique(), key=lambda x: int(x))
    cmap = plt.get_cmap("tab20")
    for i, cat in enumerate(cats):
        sub = emb[emb["cluster"] == cat]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=0.6, alpha=0.75, linewidths=0, color=cmap(i % 20), label=cat)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(f"Promoter cov10 {tag}")
    ax.legend(markerscale=5, fontsize=5, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / f"snTAPS-FL-1.promoter_cov10.coarse_{tag}.umap_by_cluster.png", dpi=600)
    plt.close(fig)

    # density图
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    hb = ax.hexbin(emb["UMAP1"], emb["UMAP2"], gridsize=90, mincnt=1)
    cb = plt.colorbar(hb, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Cell density")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(f"Cell density {tag}")
    fig.tight_layout()
    fig.savefig(outdir / f"snTAPS-FL-1.promoter_cov10.coarse_{tag}.density.png", dpi=600)
    plt.close(fig)

summary = pd.DataFrame(summary, columns=[
    "tag", "n_neighbors", "resolution", "n_cells",
    "n_clusters", "largest_cluster_n", "largest_cluster_frac", "small_cluster_frac_lt500"
])
summary.to_csv(outdir / "coarse_leiden_parameter_summary.tsv", sep="\t", index=False)
print(summary)
print("[DONE]", outdir)
