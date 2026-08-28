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

try:
    import umap
except Exception as e:
    raise SystemExit("ERROR: cannot import umap. Please install umap-learn in this env.") from e

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--global-qc", required=True)
    ap.add_argument("--bin-counts", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", default="snTAPS-FL-1")
    ap.add_argument("--min-reads", type=float, default=20000)
    ap.add_argument("--min-covered-calls", type=float, default=50000)
    ap.add_argument("--min-covered-bins", type=float, default=2000)
    ap.add_argument("--cellbin-min-cov", type=int, default=5)
    ap.add_argument("--max-cells", type=int, default=50000)
    ap.add_argument("--min-bin-cells", type=int, default=300)
    ap.add_argument("--top-variable-bins", type=int, default=2000)
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--n-clusters", type=int, default=8)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[INFO] read global qc")
    qc = pd.read_csv(args.global_qc, sep="\t")
    qc = qc[qc["CB"].astype(str) != "CB"].copy()

    for c in ["reads_used", "modified_calls", "covered_calls", "modified_rate", "covered_bins"]:
        qc[c] = pd.to_numeric(qc[c], errors="coerce")
    qc = qc.dropna()

    # high-load doublet-like upper filters
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
    cells = cells.sort_values("covered_calls", ascending=False).head(args.max_cells).copy()
    selected = set(cells["CB"])

    print(f"[INFO] total cells={len(qc)} selected cells={len(cells)}")
    cells.to_csv(outdir / f"{args.prefix}.selected_cells_for_umap.tsv", sep="\t", index=False)

    # First pass: bin statistics among selected cells
    print("[INFO] first pass: bin stats")
    bin_count = {}
    bin_sum = {}
    bin_sumsq = {}

    usecols = ["CB", "bin_id", "covered_calls", "modified_rate"]
    for chunk in pd.read_csv(args.bin_counts, sep="\t", compression="gzip", chunksize=2_000_000, usecols=usecols):
        chunk = chunk[chunk["CB"].isin(selected)]
        chunk = chunk[pd.to_numeric(chunk["covered_calls"], errors="coerce") >= args.cellbin_min_cov]
        if chunk.empty:
            continue
        chunk["modified_rate"] = pd.to_numeric(chunk["modified_rate"], errors="coerce")
        chunk = chunk.dropna(subset=["modified_rate"])
        g = chunk.groupby("bin_id")["modified_rate"].agg(["count", "sum"])
        g2 = chunk.assign(rate2=chunk["modified_rate"] ** 2).groupby("bin_id")["rate2"].sum()

        for b, row in g.iterrows():
            bin_count[b] = bin_count.get(b, 0) + int(row["count"])
            bin_sum[b] = bin_sum.get(b, 0.0) + float(row["sum"])
        for b, v in g2.items():
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

    print(f"[INFO] selected bins={len(bins)}")
    X = np.full((len(cell_list), len(bins)), np.nan, dtype=np.float32)

    print("[INFO] second pass: build matrix")
    for chunk in pd.read_csv(args.bin_counts, sep="\t", compression="gzip", chunksize=2_000_000, usecols=usecols):
        chunk = chunk[chunk["CB"].isin(selected)]
        chunk = chunk[chunk["bin_id"].isin(bin_to_j)]
        chunk = chunk[pd.to_numeric(chunk["covered_calls"], errors="coerce") >= args.cellbin_min_cov]
        if chunk.empty:
            continue
        chunk["modified_rate"] = pd.to_numeric(chunk["modified_rate"], errors="coerce")
        chunk = chunk.dropna(subset=["modified_rate"])

        for cb, b, rate in zip(chunk["CB"], chunk["bin_id"], chunk["modified_rate"]):
            X[cell_to_i[cb], bin_to_j[b]] = rate

    # impute missing with bin mean
    bin_means = np.nanmean(X, axis=0)
    inds = np.where(np.isnan(X))
    X[inds] = np.take(bin_means, inds[1])

    np.save(outdir / f"{args.prefix}.matrix_cells_by_bins.npy", X)

    print("[INFO] PCA/UMAP")
    Xs = StandardScaler(with_mean=True, with_std=True).fit_transform(X)
    pca = PCA(n_components=args.n_pcs, random_state=1)
    pcs = pca.fit_transform(Xs)

    reducer = umap.UMAP(
        n_neighbors=30,
        min_dist=0.3,
        n_components=2,
        metric="euclidean",
        random_state=1
    )
    emb = reducer.fit_transform(pcs)

    km = KMeans(n_clusters=args.n_clusters, random_state=1, n_init=20)
    cluster = km.fit_predict(pcs)

    res = cells.copy()
    res["UMAP1"] = emb[:, 0]
    res["UMAP2"] = emb[:, 1]
    res["cluster"] = cluster
    res.to_csv(outdir / f"{args.prefix}.umap_firstpass.tsv", sep="\t", index=False)

    def scatter_color(col, cmap="viridis"):
        plt.figure(figsize=(7,6))
        plt.scatter(res["UMAP1"], res["UMAP2"], c=res[col], s=3, alpha=0.7, cmap=cmap)
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

    print("[INFO] done")

if __name__ == "__main__":
    main()
