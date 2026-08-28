#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.metrics.pairwise import nan_euclidean_distances
from sklearn.cluster import KMeans
import umap


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region-counts", required=True)
    ap.add_argument("--cell-qc", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--min-cov", type=int, default=1)
    ap.add_argument("--n-top-regions", type=int, default=500)
    ap.add_argument("--min-region-cells", type=int, default=30)
    ap.add_argument("--min-regions-per-cell", type=int, default=20)
    ap.add_argument("--n-clusters", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--chunksize", type=int, default=1000000)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[INFO] read cell QC", flush=True)
    qc = pd.read_csv(args.cell_qc, sep="\t")
    for c in ["reads_used", "covered_calls", "covered_bins", "modified_rate"]:
        if c in qc.columns:
            qc[c] = pd.to_numeric(qc[c], errors="coerce")

    qc = qc.dropna(subset=["CB", "modified_rate"]).copy()
    cell_order = qc["CB"].astype(str).tolist()
    cell_set = set(cell_order)
    cell_rate = dict(zip(qc["CB"].astype(str), qc["modified_rate"]))

    print("[INFO] input cells:", len(cell_order), flush=True)

    # pass1: 统计每个promoter在这些细胞中的真实覆盖情况和残差方差
    print("[INFO] pass1 region observed stats", flush=True)
    region_count = {}
    region_sum = {}
    region_sumsq = {}
    cell_obs_all = dict.fromkeys(cell_order, 0)

    usecols = ["CB", "region_id", "modified_calls", "covered_calls", "modified_rate"]

    for chunk in pd.read_csv(args.region_counts, sep="\t", usecols=usecols,
                             chunksize=args.chunksize):
        chunk["CB"] = chunk["CB"].astype(str)
        chunk = chunk[chunk["CB"].isin(cell_set)].copy()
        if chunk.empty:
            continue

        chunk["covered_calls"] = pd.to_numeric(chunk["covered_calls"], errors="coerce")
        chunk["modified_rate"] = pd.to_numeric(chunk["modified_rate"], errors="coerce")
        chunk = chunk[(chunk["covered_calls"] >= args.min_cov) & chunk["modified_rate"].notna()].copy()
        if chunk.empty:
            continue

        # residual：区域甲基化率 - 该细胞全局甲基化率
        chunk["cell_rate"] = chunk["CB"].map(cell_rate)
        chunk["value"] = chunk["modified_rate"] - chunk["cell_rate"]

        # 每个细胞真实覆盖promoter数
        cc = chunk.groupby("CB").size()
        for cb, n in cc.items():
            cell_obs_all[cb] = cell_obs_all.get(cb, 0) + int(n)

        g = chunk.groupby("region_id")["value"].agg(["count", "sum"])
        g2 = chunk.assign(value2=chunk["value"] ** 2).groupby("region_id")["value2"].sum()

        for rid, row in g.iterrows():
            region_count[rid] = region_count.get(rid, 0) + int(row["count"])
            region_sum[rid] = region_sum.get(rid, 0.0) + float(row["sum"])
        for rid, ss in g2.items():
            region_sumsq[rid] = region_sumsq.get(rid, 0.0) + float(ss)

    stats = []
    for rid, cnt in region_count.items():
        if cnt <= 1:
            var = 0.0
        else:
            mean = region_sum[rid] / cnt
            var = region_sumsq[rid] / cnt - mean * mean
        stats.append((rid, cnt, var))

    stats = pd.DataFrame(stats, columns=["region_id", "n_cells_observed", "residual_variance"])
    stats = stats[stats["n_cells_observed"] >= args.min_region_cells].copy()
    stats = stats.sort_values(["residual_variance", "n_cells_observed"], ascending=False)
    stats = stats.head(args.n_top_regions).copy()

    stats.to_csv(outdir / f"{args.prefix}.selected_regions.tsv", sep="\t", index=False)

    selected_regions = stats["region_id"].tolist()
    selected_set = set(selected_regions)

    print("[INFO] selected regions:", len(selected_regions), flush=True)
    if len(selected_regions) < 20:
        raise ValueError("Too few observed promoter regions selected")

    # pass2: 构建含NaN的矩阵，但不做均值填补
    print("[INFO] pass2 build observed matrix with NaN", flush=True)
    cell_index = {cb: i for i, cb in enumerate(cell_order)}
    reg_index = {r: j for j, r in enumerate(selected_regions)}
    X = np.full((len(cell_order), len(selected_regions)), np.nan, dtype=np.float32)

    for chunk in pd.read_csv(args.region_counts, sep="\t", usecols=usecols,
                             chunksize=args.chunksize):
        chunk["CB"] = chunk["CB"].astype(str)
        chunk = chunk[chunk["CB"].isin(cell_set) & chunk["region_id"].isin(selected_set)].copy()
        if chunk.empty:
            continue

        chunk["covered_calls"] = pd.to_numeric(chunk["covered_calls"], errors="coerce")
        chunk["modified_rate"] = pd.to_numeric(chunk["modified_rate"], errors="coerce")
        chunk = chunk[(chunk["covered_calls"] >= args.min_cov) & chunk["modified_rate"].notna()].copy()
        if chunk.empty:
            continue

        chunk["cell_rate"] = chunk["CB"].map(cell_rate)
        chunk["value"] = chunk["modified_rate"] - chunk["cell_rate"]

        # 如果同一个CB-region有重复，取均值
        chunk = chunk.groupby(["CB", "region_id"], as_index=False)["value"].mean()

        for cb, rid, val in chunk[["CB", "region_id", "value"]].itertuples(index=False):
            X[cell_index[cb], reg_index[rid]] = val

    obs_per_cell = np.isfinite(X).sum(axis=1)
    keep = obs_per_cell >= args.min_regions_per_cell
    X = X[keep, :]
    qc2 = qc.loc[keep].copy()
    obs_per_cell = obs_per_cell[keep]

    print("[INFO] cells after observed promoter filter:", X.shape[0], flush=True)
    print("[INFO] matrix shape:", X.shape, flush=True)

    nan_fraction = np.isnan(X).mean()
    print("[INFO] nan_fraction_no_fill:", nan_fraction, flush=True)

    if X.shape[0] < 50:
        raise ValueError("Too few cells after promoter observed filter")

    # nan_euclidean_distances：只用共同非NaN位置计算距离
    print("[INFO] compute nan-aware distances", flush=True)
    D = nan_euclidean_distances(X)

    finite = np.isfinite(D)
    maxd = np.nanmax(D[finite]) if finite.any() else 1.0
    D[~finite] = maxd * 1.2
    np.fill_diagonal(D, 0.0)

    print("[INFO] UMAP from precomputed observed-only distance", flush=True)
    reducer = umap.UMAP(
        n_neighbors=min(30, X.shape[0] - 1),
        min_dist=0.3,
        metric="precomputed",
        random_state=args.seed
    )
    emb = reducer.fit_transform(D)

    n_clusters = min(args.n_clusters, max(2, X.shape[0] // 50))
    cluster = KMeans(n_clusters=n_clusters, random_state=args.seed, n_init=20).fit_predict(emb)

    res = qc2.copy()
    res["UMAP1"] = emb[:, 0]
    res["UMAP2"] = emb[:, 1]
    res["cluster"] = cluster
    res["observed_selected_promoters"] = obs_per_cell

    res.to_csv(outdir / f"{args.prefix}.umap.tsv", sep="\t", index=False)
    res["cluster"].value_counts().sort_index().to_csv(
        outdir / f"{args.prefix}.cluster_counts.tsv",
        sep="\t",
        header=["n_cells"]
    )

    matrix_qc = pd.DataFrame([{
        "n_cells_input": len(cell_order),
        "n_cells_used": X.shape[0],
        "n_regions_selected": X.shape[1],
        "min_cov": args.min_cov,
        "min_region_cells": args.min_region_cells,
        "min_regions_per_cell": args.min_regions_per_cell,
        "nan_fraction_no_fill": nan_fraction,
        "observed_promoters_median": float(np.median(obs_per_cell)),
        "observed_promoters_q10": float(np.quantile(obs_per_cell, 0.10)),
        "observed_promoters_q90": float(np.quantile(obs_per_cell, 0.90)),
    }])
    matrix_qc.to_csv(outdir / f"{args.prefix}.matrix_qc.tsv", sep="\t", index=False)

    plt.figure(figsize=(7, 6))
    sc = plt.scatter(res["UMAP1"], res["UMAP2"], c=res["cluster"], s=4, cmap="tab10")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(f"{args.prefix} observed-only promoter")
    plt.colorbar(sc, label="cluster")
    plt.tight_layout()
    plt.savefig(outdir / f"{args.prefix}.umap_by_cluster.png", dpi=220)
    plt.close()

    for c in ["reads_used", "covered_calls", "covered_bins", "modified_rate", "observed_selected_promoters"]:
        if c not in res.columns:
            continue
        plt.figure(figsize=(7, 6))
        sc = plt.scatter(res["UMAP1"], res["UMAP2"], c=pd.to_numeric(res[c], errors="coerce"), s=4)
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
