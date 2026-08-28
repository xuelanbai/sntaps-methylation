#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import umap

def logit(x):
    x = np.clip(x, 1e-4, 1 - 1e-4)
    return np.log(x / (1 - x))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--global-qc", required=True)
    ap.add_argument("--bin-counts", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--mode", choices=["residual", "seurat_norm", "raw"], required=True)
    ap.add_argument("--feature-bin-list", default=None)
    ap.add_argument("--min-reads", type=float, default=20000)
    ap.add_argument("--min-covered-calls", type=float, default=50000)
    ap.add_argument("--min-covered-bins", type=float, default=2000)
    ap.add_argument("--cellbin-min-cov", type=int, default=5)
    ap.add_argument("--max-cells", type=int, default=50000)
    ap.add_argument("--min-bin-cells", type=int, default=300)
    ap.add_argument("--top-variable-bins", type=int, default=2000)
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--n-clusters", type=int, default=8)
    ap.add_argument("--n-jobs", type=int, default=16)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    feature_bins = None
    if args.feature_bin_list:
        with open(args.feature_bin_list) as f:
            feature_bins = set(x.strip().split()[0] for x in f if x.strip())
        print(f"[INFO] feature bins loaded: {len(feature_bins)}", flush=True)

    print("[INFO] read global qc", flush=True)
    qc = pd.read_csv(args.global_qc, sep="\t")
    qc = qc[qc["CB"].astype(str) != "CB"].copy()

    for c in ["reads_used", "modified_calls", "covered_calls", "modified_rate", "covered_bins"]:
        qc[c] = pd.to_numeric(qc[c], errors="coerce")
    qc = qc.dropna().copy()

    reads_p995 = qc["reads_used"].quantile(0.995)
    covcalls_p995 = qc["covered_calls"].quantile(0.995)
    covbins_p995 = qc["covered_bins"].quantile(0.995)

    keep = (
        (qc["reads_used"] >= args.min_reads) &
        (qc["covered_calls"] >= args.min_covered_calls) &
        (qc["covered_bins"] >= args.min_covered_bins) &
        (qc["reads_used"] <= reads_p995) &
        (qc["covered_calls"] <= covcalls_p995) &
        (qc["covered_bins"] <= covbins_p995)
    )

    cells = qc.loc[keep].copy()
    if len(cells) > args.max_cells:
        cells = cells.sample(n=args.max_cells, random_state=args.seed).copy()
    cells = cells.reset_index(drop=True)

    selected = set(cells["CB"])
    global_rate = dict(zip(cells["CB"], cells["modified_rate"]))

    print(f"[INFO] total cells={len(qc)} qc_pass={keep.sum()} selected={len(cells)}", flush=True)
    cells.to_csv(outdir / f"{args.prefix}.selected_cells.tsv", sep="\t", index=False)

    usecols = ["CB", "bin_id", "modified_calls", "covered_calls"]

    bin_count, bin_sum, bin_sumsq = {}, {}, {}

    print("[INFO] first pass: bin variance", flush=True)
    for chunk in pd.read_csv(args.bin_counts, sep="\t", compression="gzip", chunksize=2000000, usecols=usecols):
        chunk = chunk[chunk["CB"].isin(selected)].copy()
        if feature_bins is not None:
            chunk = chunk[chunk["bin_id"].isin(feature_bins)].copy()
        if chunk.empty:
            continue

        chunk["modified_calls"] = pd.to_numeric(chunk["modified_calls"], errors="coerce")
        chunk["covered_calls"] = pd.to_numeric(chunk["covered_calls"], errors="coerce")
        chunk = chunk.dropna()
        chunk = chunk[chunk["covered_calls"] >= args.cellbin_min_cov].copy()
        if chunk.empty:
            continue

        m = chunk["modified_calls"].to_numpy(dtype=float)
        c = chunk["covered_calls"].to_numpy(dtype=float)
        rate = (m + 0.5) / (c + 1.0)

        if args.mode == "residual":
            g = chunk["CB"].map(global_rate).astype(float).to_numpy()
            val = logit(rate) - logit(g)
            val = np.clip(val, -5, 5)
        else:
            val = rate

        chunk["value"] = val

        stat = chunk.groupby("bin_id")["value"].agg(["count", "sum"])
        stat2 = chunk.assign(v2=chunk["value"] ** 2).groupby("bin_id")["v2"].sum()

        for b, row in stat.iterrows():
            bin_count[b] = bin_count.get(b, 0) + int(row["count"])
            bin_sum[b] = bin_sum.get(b, 0.0) + float(row["sum"])
        for b, v in stat2.items():
            bin_sumsq[b] = bin_sumsq.get(b, 0.0) + float(v)

    rows = []
    for b, cnt in bin_count.items():
        if cnt >= args.min_bin_cells:
            mean = bin_sum[b] / cnt
            var = max(bin_sumsq[b] / cnt - mean * mean, 0.0)
            rows.append((b, cnt, mean, var))

    bin_stats = pd.DataFrame(rows, columns=["bin_id", "n_cells", "mean", "var"])
    bin_stats = bin_stats.sort_values("var", ascending=False).head(args.top_variable_bins).copy()
    bin_stats.to_csv(outdir / f"{args.prefix}.selected_variable_bins.tsv", sep="\t", index=False)

    bins = list(bin_stats["bin_id"])
    bin_to_j = {b: j for j, b in enumerate(bins)}
    cell_list = list(cells["CB"])
    cell_to_i = {c: i for i, c in enumerate(cell_list)}

    print(f"[INFO] selected bins={len(bins)}", flush=True)

    fill_value = 0.0 if args.mode == "residual" else np.nan
    X = np.full((len(cell_list), len(bins)), fill_value, dtype=np.float32)

    print("[INFO] second pass: build matrix", flush=True)
    for chunk in pd.read_csv(args.bin_counts, sep="\t", compression="gzip", chunksize=2000000, usecols=usecols):
        chunk = chunk[chunk["CB"].isin(selected)].copy()
        chunk = chunk[chunk["bin_id"].isin(bin_to_j)].copy()
        if chunk.empty:
            continue

        chunk["modified_calls"] = pd.to_numeric(chunk["modified_calls"], errors="coerce")
        chunk["covered_calls"] = pd.to_numeric(chunk["covered_calls"], errors="coerce")
        chunk = chunk.dropna()
        chunk = chunk[chunk["covered_calls"] >= args.cellbin_min_cov].copy()
        if chunk.empty:
            continue

        m = chunk["modified_calls"].to_numpy(dtype=float)
        c = chunk["covered_calls"].to_numpy(dtype=float)
        rate = (m + 0.5) / (c + 1.0)

        if args.mode == "residual":
            g = chunk["CB"].map(global_rate).astype(float).to_numpy()
            val = logit(rate) - logit(g)
            val = np.clip(val, -5, 5)
        else:
            val = rate

        for cb, b, v in zip(chunk["CB"], chunk["bin_id"], val):
            X[cell_to_i[cb], bin_to_j[b]] = v

    if args.mode in ["raw", "seurat_norm"]:
        bin_means = np.nanmean(X, axis=0)
        inds = np.where(np.isnan(X))
        X[inds] = np.take(bin_means, inds[1])

    if args.mode == "seurat_norm":
        row_sum = X.sum(axis=1)
        row_sum[row_sum == 0] = 1
        X = np.log1p(X / row_sum[:, None] * 10000).astype(np.float32)

    print("[INFO] PCA / regress QC covariates / UMAP", flush=True)
    Xs = StandardScaler(with_mean=True, with_std=True).fit_transform(X)
    pcs = PCA(n_components=args.n_pcs, random_state=args.seed).fit_transform(Xs)

    qc_mat = np.column_stack([
        np.log10(cells["reads_used"].to_numpy(dtype=float) + 1.0),
        np.log10(cells["covered_calls"].to_numpy(dtype=float) + 1.0),
        np.log10(cells["covered_bins"].to_numpy(dtype=float) + 1.0),
        cells["modified_rate"].to_numpy(dtype=float)
    ])

    qc_mat = StandardScaler(with_mean=True, with_std=True).fit_transform(qc_mat)
    design = np.column_stack([np.ones(qc_mat.shape[0]), qc_mat])

    beta = np.linalg.lstsq(design, pcs, rcond=None)[0]
    fitted_qc = qc_mat @ beta[1:, :]
    pcs_for_umap = pcs - fitted_qc

    reducer = umap.UMAP(
        n_neighbors=30,
        min_dist=0.3,
        n_components=2,
        metric="euclidean",
        random_state=None,
        n_jobs=args.n_jobs
    )
    emb = reducer.fit_transform(pcs_for_umap)

    cluster = KMeans(n_clusters=args.n_clusters, random_state=args.seed, n_init=20).fit_predict(pcs_for_umap)

    res = cells.copy()
    res["UMAP1"] = emb[:, 0]
    res["UMAP2"] = emb[:, 1]
    res["cluster"] = cluster
    res.to_csv(outdir / f"{args.prefix}.umap.tsv", sep="\t", index=False)

    corr = res[["UMAP1", "UMAP2", "reads_used", "covered_calls", "covered_bins", "modified_rate"]].corr(method="spearman")
    corr.to_csv(outdir / f"{args.prefix}.umap_qc_spearman_corr.tsv", sep="\t")

    pc_qc = pd.DataFrame(pcs_for_umap[:, :min(10, pcs_for_umap.shape[1])],
                         columns=[f"regressed_PC{i+1}" for i in range(min(10, pcs_for_umap.shape[1]))])
    pc_qc["reads_used"] = cells["reads_used"].to_numpy()
    pc_qc["covered_calls"] = cells["covered_calls"].to_numpy()
    pc_qc["covered_bins"] = cells["covered_bins"].to_numpy()
    pc_qc["modified_rate"] = cells["modified_rate"].to_numpy()
    pc_corr = pc_qc.corr(method="spearman")
    pc_corr.to_csv(outdir / f"{args.prefix}.regressed_PC_qc_spearman_corr.tsv", sep="\t")

    def scatter_color(col):
        plt.figure(figsize=(7,6))
        plt.scatter(res["UMAP1"], res["UMAP2"], c=res[col], s=3, alpha=0.7)
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.title(f"{args.prefix} colored by {col}")
        cb = plt.colorbar()
        cb.set_label(col)
        plt.tight_layout()
        plt.savefig(outdir / f"{args.prefix}.umap_by_{col}.png", dpi=200)
        plt.close()

    for col in ["reads_used", "covered_calls", "covered_bins", "modified_rate"]:
        scatter_color(col)

    plt.figure(figsize=(7,6))
    plt.scatter(res["UMAP1"], res["UMAP2"], c=res["cluster"], s=3, alpha=0.7, cmap="tab10")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(f"{args.prefix} KMeans clusters")
    cb = plt.colorbar()
    cb.set_label("cluster")
    plt.tight_layout()
    plt.savefig(outdir / f"{args.prefix}.umap_by_cluster.png", dpi=200)
    plt.close()

    print("[INFO] done", flush=True)

if __name__ == "__main__":
    main()
