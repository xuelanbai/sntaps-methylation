#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin-counts", required=True)
    ap.add_argument("--global-cb", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--min-cell-covered-cpg", type=int, default=5000)
    ap.add_argument("--min-cell-covered-bins", type=int, default=200)
    ap.add_argument("--min-cellbin-cov", type=int, default=5)
    ap.add_argument("--min-bin-cells", type=int, default=20)
    ap.add_argument("--top-variable-bins", type=int, default=2000)
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--n-clusters", type=int, default=8)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[INFO] loading global CB table")
    g = pd.read_csv(args.global_cb, sep="\t")

    g_filt = g[
        (g["covered_calls"] >= args.min_cell_covered_cpg) &
        (g["covered_bins"] >= args.min_cell_covered_bins)
    ].copy()

    keep_cells = set(g_filt["CB"])
    print(f"[INFO] cells before filter: {g.shape[0]}")
    print(f"[INFO] cells after global filter: {len(keep_cells)}")

    if len(keep_cells) < 10:
        raise SystemExit("[ERROR] Too few cells after global filter. Lower min-cell thresholds.")

    print("[INFO] loading bin count table")
    d = pd.read_csv(args.bin_counts, sep="\t")

    d = d[d["CB"].isin(keep_cells)]
    d = d[d["covered_calls"] >= args.min_cellbin_cov]

    print(f"[INFO] cell-bin rows after covered_calls filter: {d.shape[0]}")

    bin_cell_n = d.groupby("bin_id")["CB"].nunique()
    keep_bins = set(bin_cell_n[bin_cell_n >= args.min_bin_cells].index)
    d = d[d["bin_id"].isin(keep_bins)]

    print(f"[INFO] bins after min-bin-cells filter: {len(keep_bins)}")
    print(f"[INFO] rows after bin filter: {d.shape[0]}")

    if d.shape[0] == 0:
        raise SystemExit("[ERROR] No cell-bin rows left. Lower min-cellbin-cov or min-bin-cells.")

    print("[INFO] making cell x bin matrix")
    mat = d.pivot_table(
        index="CB",
        columns="bin_id",
        values="modified_rate",
        aggfunc="mean"
    )

    print(f"[INFO] matrix shape before imputation: {mat.shape}")

    if mat.shape[0] < 10 or mat.shape[1] < 20:
        raise SystemExit("[ERROR] Matrix too small for clustering. Lower thresholds or use more reads.")

    # bin mean imputation
    mat = mat.apply(lambda x: x.fillna(x.mean()), axis=0)
    mat = mat.fillna(mat.stack().mean())

    # variable bins
    variances = mat.var(axis=0).sort_values(ascending=False)
    topn = min(args.top_variable_bins, len(variances))
    top_bins = list(variances.head(topn).index)
    mat = mat[top_bins]

    print(f"[INFO] matrix shape after variable-bin selection: {mat.shape}")

    mat.to_csv(outdir / f"{args.prefix}.cell_by_bin_matrix.tsv", sep="\t")
    variances.to_csv(outdir / f"{args.prefix}.bin_variance.tsv", sep="\t", header=["variance"])

    # scale features
    X = mat.values.astype(float)
    X = X - X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, keepdims=True)
    sd[sd == 0] = 1
    X = X / sd
    X = np.clip(X, -10, 10)

    from sklearn.decomposition import PCA
    from sklearn.cluster import KMeans
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_pcs = min(args.n_pcs, X.shape[0] - 1, X.shape[1] - 1)
    if n_pcs < 2:
        raise SystemExit("[ERROR] Too few cells/features for PCA.")

    print(f"[INFO] PCA n_components={n_pcs}")
    pca = PCA(n_components=n_pcs, random_state=1)
    X_pca = pca.fit_transform(X)

    pd.DataFrame(
        X_pca,
        index=mat.index,
        columns=[f"PC{i+1}" for i in range(n_pcs)]
    ).to_csv(outdir / f"{args.prefix}.pca.tsv", sep="\t")

    pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(n_pcs)],
        "explained_variance_ratio": pca.explained_variance_ratio_
    }).to_csv(outdir / f"{args.prefix}.pca_variance.tsv", sep="\t", index=False)

    n_clusters = min(args.n_clusters, max(2, X_pca.shape[0] // 20))
    n_clusters = max(2, n_clusters)

    print(f"[INFO] KMeans n_clusters={n_clusters}")
    km = KMeans(n_clusters=n_clusters, random_state=1, n_init=20)
    clusters = km.fit_predict(X_pca)

    # UMAP if available
    use_umap = True
    try:
        import umap
        reducer = umap.UMAP(
            n_neighbors=min(15, X_pca.shape[0] - 1),
            min_dist=0.3,
            metric="euclidean",
            random_state=1
        )
        emb = reducer.fit_transform(X_pca)
        dim1_name, dim2_name = "UMAP1", "UMAP2"
        plot_name = "umap"
        print("[INFO] UMAP finished")
    except Exception as e:
        print("[WARN] UMAP unavailable or failed; using PCA1/PCA2 instead.")
        print("[WARN]", str(e))
        emb = X_pca[:, :2]
        dim1_name, dim2_name = "PC1", "PC2"
        plot_name = "pca"

    obs = pd.DataFrame({
        "CB": mat.index,
        dim1_name: emb[:, 0],
        dim2_name: emb[:, 1],
        "cluster": clusters.astype(str)
    })

    obs = obs.merge(
        g[["CB", "reads_used", "modified_calls", "covered_calls", "modified_rate", "covered_bins"]],
        on="CB",
        how="left"
    )

    obs.to_csv(outdir / f"{args.prefix}.{plot_name}_clusters.tsv", sep="\t", index=False)

    cluster_summary = obs.groupby("cluster").agg(
        n_cells=("CB", "count"),
        median_reads=("reads_used", "median"),
        median_covered_calls=("covered_calls", "median"),
        median_modified_rate=("modified_rate", "median"),
        median_covered_bins=("covered_bins", "median")
    ).reset_index()

    cluster_summary.to_csv(outdir / f"{args.prefix}.cluster_summary.tsv", sep="\t", index=False)

    plt.figure(figsize=(6, 5))
    for cl in sorted(obs["cluster"].unique(), key=lambda x: int(x)):
        sub = obs[obs["cluster"] == cl]
        plt.scatter(sub[dim1_name], sub[dim2_name], s=6, alpha=0.8, label=cl)

    plt.xlabel(dim1_name)
    plt.ylabel(dim2_name)
    plt.title(args.prefix)
    plt.legend(title="cluster", markerscale=2, fontsize=7)
    plt.tight_layout()
    plt.savefig(outdir / f"{args.prefix}.{plot_name}_clusters.png", dpi=300)
    plt.close()

    print(f"[INFO] wrote {outdir / f'{args.prefix}.{plot_name}_clusters.png'}")
    print("[INFO] done")

if __name__ == "__main__":
    main()
