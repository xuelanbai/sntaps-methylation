#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import umap

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region-counts", required=True)
    ap.add_argument("--cell-qc", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)

    ap.add_argument("--min-region-cells", type=int, default=100)
    ap.add_argument("--n-top-regions", type=int, default=3000)
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--n-clusters", type=int, default=4)
    ap.add_argument("--n-neighbors", type=int, default=30)
    ap.add_argument("--min-dist", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--qc-regress", default="reads_used,modified_rate")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[INFO] read cell QC", flush=True)
    qc = pd.read_csv(args.cell_qc, sep="\t")
    qc["CB"] = qc["CB"].astype(str)

    for c in ["reads_used", "covered_calls", "covered_bins", "modified_rate"]:
        if c in qc.columns:
            qc[c] = pd.to_numeric(qc[c], errors="coerce")

    keep_cells = set(qc["CB"])
    print("[INFO] n cells in QC:", len(keep_cells), flush=True)

    print("[INFO] read region counts", flush=True)
    df = pd.read_csv(args.region_counts, sep="\t")
    df["CB"] = df["CB"].astype(str)
    df = df[df["CB"].isin(keep_cells)].copy()

    region_col = "region_id" if "region_id" in df.columns else df.columns[1]
    df[region_col] = df[region_col].astype(str)

    df["modified_calls"] = pd.to_numeric(df["modified_calls"], errors="coerce").fillna(0)
    df["covered_calls"] = pd.to_numeric(df["covered_calls"], errors="coerce").fillna(0)
    df = df[df["covered_calls"] > 0].copy()

    print("[INFO] observed rows:", len(df), flush=True)

    print("[INFO] region observed stats", flush=True)
    stats = df.groupby(region_col).agg(
        n_cells=("CB", "nunique"),
        total_cov=("covered_calls", "sum"),
        total_mod=("modified_calls", "sum"),
    )
    stats["region_rate"] = stats["total_mod"] / stats["total_cov"].replace(0, np.nan)
    stats = stats.replace([np.inf, -np.inf], np.nan).fillna(0)

    stats.to_csv(outdir / f"{args.prefix}.all_region_observed_stats.tsv", sep="\t")

    selected = stats[stats["n_cells"] >= args.min_region_cells].copy()
    selected = selected.sort_values(["n_cells", "total_cov"], ascending=False)
    selected = selected.head(args.n_top_regions)

    print("[INFO] selected regions:", len(selected), flush=True)
    if len(selected) < 10:
        raise ValueError(f"too few selected regions: {len(selected)}")

    selected.to_csv(outdir / f"{args.prefix}.selected_region_stats.tsv", sep="\t")

    sel_regions = list(selected.index)
    region_to_i = {r:i for i,r in enumerate(sel_regions)}
    cell_to_i = {c:i for i,c in enumerate(qc["CB"])}

    df = df[df[region_col].isin(region_to_i)].copy()
    print("[INFO] rows after selected regions:", len(df), flush=True)

    cell_rate = qc.set_index("CB")["modified_rate"].to_dict()
    global_rate = float(qc["modified_rate"].median())

    rows = df["CB"].map(cell_to_i).to_numpy()
    cols = df[region_col].map(region_to_i).to_numpy()

    cov = df["covered_calls"].to_numpy(dtype=np.float64)
    mod = df["modified_calls"].to_numpy(dtype=np.float64)
    rate = df["CB"].map(cell_rate).fillna(global_rate).to_numpy(dtype=np.float64)
    rate = np.clip(rate, 1e-4, 0.9999)

    # coverage-weighted residual: observed - expected, scaled by binomial variance
    expected = cov * rate
    denom = np.sqrt(cov * rate * (1 - rate) + 1.0)
    val = (mod - expected) / denom
    val = np.clip(val, -10, 10).astype(np.float32)

    X = sparse.csr_matrix((val, (rows, cols)), shape=(len(qc), len(sel_regions)))
    obs_fraction = X.nnz / (X.shape[0] * X.shape[1])

    print("[INFO] matrix shape:", X.shape, flush=True)
    print("[INFO] observed fraction:", obs_fraction, flush=True)

    n_components = min(args.n_pcs, X.shape[0] - 1, X.shape[1] - 1)
    if n_components < 2:
        raise ValueError(f"too few PCs possible: {n_components}")

    print("[INFO] scale sparse matrix", flush=True)
    Xs = StandardScaler(with_mean=False).fit_transform(X)

    print("[INFO] TruncatedSVD", flush=True)
    pcs = TruncatedSVD(n_components=n_components, random_state=args.seed).fit_transform(Xs)

    # 回归QC，只在PC空间做，避免直接破坏稀疏矩阵
    qc_keys = [x.strip() for x in args.qc_regress.split(",") if x.strip()]
    qc_keys = [x for x in qc_keys if x in qc.columns]

    if qc_keys:
        print("[INFO] regress QC from PCs:", qc_keys, flush=True)
        Z = qc[qc_keys].copy()
        for c in qc_keys:
            Z[c] = pd.to_numeric(Z[c], errors="coerce")
            Z[c] = Z[c].fillna(Z[c].median())
            if Z[c].max() > 10:
                Z[c] = np.log10(Z[c] + 1)

        Z = StandardScaler().fit_transform(Z)
        design = np.column_stack([np.ones(Z.shape[0]), Z])
        beta = np.linalg.lstsq(design, pcs, rcond=None)[0]
        pcs = pcs - design[:, 1:] @ beta[1:, :]

    print("[INFO] UMAP", flush=True)
    nn = min(args.n_neighbors, len(qc) - 1)
    emb = umap.UMAP(
        n_neighbors=nn,
        min_dist=args.min_dist,
        metric="euclidean",
        random_state=args.seed,
        n_components=2,
    ).fit_transform(pcs)

    print("[INFO] KMeans", flush=True)
    k = min(args.n_clusters, len(qc))
    cluster = KMeans(n_clusters=k, random_state=args.seed, n_init=20).fit_predict(pcs)

    res = qc.copy()
    res["UMAP1"] = emb[:, 0]
    res["UMAP2"] = emb[:, 1]
    res["cluster"] = cluster

    res.to_csv(outdir / f"{args.prefix}.umap.tsv", sep="\t", index=False)
    res["cluster"].value_counts().sort_index().to_csv(
        outdir / f"{args.prefix}.cluster_counts.tsv",
        sep="\t",
        header=["n_cells"]
    )

    qc_plot_cols = [c for c in ["reads_used", "covered_calls", "covered_bins", "modified_rate"] if c in res.columns]
    corr_cols = ["UMAP1", "UMAP2"] + qc_plot_cols
    res[corr_cols].corr(method="spearman").to_csv(
        outdir / f"{args.prefix}.umap_qc_spearman_corr.tsv",
        sep="\t"
    )

    pd.DataFrame([{
        "n_cells": X.shape[0],
        "n_regions": X.shape[1],
        "observed_rows": int(X.nnz),
        "observed_fraction": obs_fraction,
        "min_region_cells": args.min_region_cells,
        "n_top_regions": args.n_top_regions,
        "n_pcs_used": n_components,
        "qc_regress": ",".join(qc_keys),
    }]).to_csv(outdir / f"{args.prefix}.matrix_qc.tsv", sep="\t", index=False)

    plt.figure(figsize=(7,6))
    sc = plt.scatter(res["UMAP1"], res["UMAP2"], c=res["cluster"], s=4, cmap="tab10")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(f"{args.prefix} sparse weighted residual")
    plt.colorbar(sc, label="cluster")
    plt.tight_layout()
    plt.savefig(outdir / f"{args.prefix}.umap_by_cluster.png", dpi=220)
    plt.close()

    for c in qc_plot_cols:
        plt.figure(figsize=(7,6))
        sc = plt.scatter(res["UMAP1"], res["UMAP2"], c=res[c], s=4)
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.title(f"{args.prefix} colored by {c}")
        plt.colorbar(sc, label=c)
        plt.tight_layout()
        plt.savefig(outdir / f"{args.prefix}.umap_by_{c}.png", dpi=220)
        plt.close()

    print("[INFO] done", flush=True)

if __name__ == "__main__":
    main()
