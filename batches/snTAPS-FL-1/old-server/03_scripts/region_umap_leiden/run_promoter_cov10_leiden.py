import pandas as pd
import numpy as np
from pathlib import Path
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import scanpy as sc
import anndata as ad

base = Path("/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation")

region_counts = base / "08_qc_summary/final_single_sample/region_matrix_covfilter/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.cov10.region_counts.tsv.gz"
old_umap_dir = base / "08_qc_summary/final_single_sample/region_umap_paper/P10_promoter_residual_cov10_dropPC1_smooth100"
selected_regions_file = old_umap_dir / "snTAPS-FL-1.promoter_residual_cov10_dropPC1_smooth100.selected_regions.tsv"
old_umap_file = old_umap_dir / "snTAPS-FL-1.promoter_residual_cov10_dropPC1_smooth100.umap.tsv"

outdir = base / "08_qc_summary/final_single_sample/region_umap_leiden/promoter_cov10_leiden"
outdir.mkdir(parents=True, exist_ok=True)

print("[INFO] reading selected cells and regions")
cells_df = pd.read_csv(old_umap_file, sep="\t", usecols=["CB", "reads_used", "covered_calls", "covered_bins", "modified_rate"])
cells = cells_df["CB"].astype(str).tolist()
cell_to_i = {c:i for i,c in enumerate(cells)}

regions_df = pd.read_csv(selected_regions_file, sep="\t")
region_col = "region_id" if "region_id" in regions_df.columns else regions_df.columns[0]
regions = regions_df[region_col].astype(str).tolist()
region_to_j = {r:j for j,r in enumerate(regions)}

print("[INFO] cells:", len(cells))
print("[INFO] selected regions:", len(regions))

rows = []
cols = []
vals = []

print("[INFO] reading region_counts in chunks")
usecols = ["CB", "region_id", "modified_rate"]
for chunk in pd.read_csv(region_counts, sep="\t", compression="gzip", usecols=usecols, dtype={"CB":str, "region_id":str}, chunksize=1000000):
    chunk = chunk[chunk["CB"].isin(cell_to_i) & chunk["region_id"].isin(region_to_j)]
    if len(chunk) == 0:
        continue
    rows.extend(chunk["CB"].map(cell_to_i).to_numpy())
    cols.extend(chunk["region_id"].map(region_to_j).to_numpy())
    vals.extend(pd.to_numeric(chunk["modified_rate"], errors="coerce").fillna(0).to_numpy())

print("[INFO] nonzero pairs:", len(vals))

X = sparse.coo_matrix(
    (np.asarray(vals, dtype=np.float32), (np.asarray(rows), np.asarray(cols))),
    shape=(len(cells), len(regions))
).tocsr()

# 缺失不是0甲基化。这里先用region均值填补缺失，再做中心化。
print("[INFO] filling missing values by region mean")
sum_by_region = np.asarray(X.sum(axis=0)).ravel()
cnt_by_region = np.asarray((X > 0).sum(axis=0)).ravel()
mean_by_region = np.divide(sum_by_region, cnt_by_region, out=np.zeros_like(sum_by_region, dtype=float), where=cnt_by_region>0)
mean_by_region = np.nan_to_num(mean_by_region, nan=0.0)

Xd = np.tile(mean_by_region, (X.shape[0], 1)).astype(np.float32)
coo = X.tocoo()
Xd[coo.row, coo.col] = coo.data

# region中心化
Xd = Xd - Xd.mean(axis=0, keepdims=True)

# 标准化，避免少数region支配
Xd = StandardScaler(with_mean=True, with_std=True).fit_transform(Xd).astype(np.float32)

print("[INFO] PCA/SVD")
svd = TruncatedSVD(n_components=50, random_state=1)
pcs = svd.fit_transform(Xd).astype(np.float32)

# 回归QC
qc_cols = ["reads_used", "covered_calls", "covered_bins", "modified_rate"]
qc = cells_df[qc_cols].copy()
for c in qc_cols:
    qc[c] = pd.to_numeric(qc[c], errors="coerce")
qc = np.log1p(qc.fillna(qc.median()).to_numpy())
qc = StandardScaler().fit_transform(qc)

print("[INFO] regress QC from PCs")
pcs_reg = pcs.copy()
lr = LinearRegression()
for k in range(pcs.shape[1]):
    lr.fit(qc, pcs[:, k])
    pcs_reg[:, k] = pcs[:, k] - lr.predict(qc)

# drop PC1
pcs_use = pcs_reg[:, 1:31]

adata = ad.AnnData(X=pcs_use)
adata.obs["CB"] = cells
adata.obs_names = cells
for c in qc_cols:
    adata.obs[c] = cells_df[c].to_numpy()

print("[INFO] neighbors/UMAP/Leiden")
sc.pp.neighbors(adata, n_neighbors=60, n_pcs=30, use_rep="X")
sc.tl.umap(adata, min_dist=0.35, random_state=1)

for res in [0.3, 0.5, 0.8, 1.0]:
    key = f"leiden_{res}"
    sc.tl.leiden(adata, resolution=res, key_added=key)

emb = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs_names)
emb.insert(0, "CB", emb.index)
for c in qc_cols:
    emb[c] = adata.obs[c].values
for res in [0.3, 0.5, 0.8, 1.0]:
    emb[f"leiden_{res}"] = adata.obs[f"leiden_{res}"].astype(str).values

emb.to_csv(outdir / "snTAPS-FL-1.promoter_cov10.leiden_umap.tsv", sep="\t", index=False)

# 相关性
corr = emb[["UMAP1", "UMAP2"] + qc_cols].corr(method="spearman")
corr.to_csv(outdir / "snTAPS-FL-1.promoter_cov10.leiden_umap_qc_spearman_corr.tsv", sep="\t")
print("[INFO] QC corr")
print(corr.loc[["UMAP1","UMAP2"], qc_cols])

# 画图
def plot_cat(col, filename):
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    cats = sorted(emb[col].unique(), key=lambda x: int(x) if str(x).isdigit() else str(x))
    cmap = plt.get_cmap("tab20")
    for i, cat in enumerate(cats):
        sub = emb[emb[col] == cat]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=0.6, alpha=0.75, linewidths=0, label=str(cat), color=cmap(i % 20))
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(col)
    ax.legend(markerscale=5, fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / filename, dpi=600)
    plt.close(fig)

def plot_cont(col, filename):
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    sca = ax.scatter(emb["UMAP1"], emb["UMAP2"], c=emb[col], s=0.6, alpha=0.75, linewidths=0)
    cb = plt.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(col)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(col)
    fig.tight_layout()
    fig.savefig(outdir / filename, dpi=600)
    plt.close(fig)

for res in [0.3, 0.5, 0.8, 1.0]:
    plot_cat(f"leiden_{res}", f"snTAPS-FL-1.promoter_cov10.umap_by_leiden_res{res}.png")

for c in qc_cols:
    plot_cont(c, f"snTAPS-FL-1.promoter_cov10.umap_by_{c}.png")

# cluster count
for res in [0.3, 0.5, 0.8, 1.0]:
    tab = emb[f"leiden_{res}"].value_counts().sort_index()
    tab.to_csv(outdir / f"snTAPS-FL-1.promoter_cov10.leiden_res{res}.cluster_counts.tsv", sep="\t")

adata.write_h5ad(outdir / "snTAPS-FL-1.promoter_cov10.leiden.h5ad")
print("[INFO] done:", outdir)
