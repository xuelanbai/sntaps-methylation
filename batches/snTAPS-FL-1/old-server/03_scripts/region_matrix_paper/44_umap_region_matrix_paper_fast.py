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
    ap.add_argument("--region-counts", required=True)
    ap.add_argument("--global-qc", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)

    ap.add_argument("--mode", choices=["residual", "logit_rate", "activity_5mc"], default="residual")
    ap.add_argument("--max-cells", type=int, default=50000)
    ap.add_argument("--min-cell-regions", type=int, default=500)
    ap.add_argument("--min-region-cells", type=int, default=500)
    ap.add_argument("--top-variable-regions", type=int, default=2000)

    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--drop-pc1", action="store_true")
    ap.add_argument("--regress-qc", action="store_true")
    ap.add_argument("--n-clusters", type=int, default=8)
    ap.add_argument("--n-neighbors", type=int, default=30)
    ap.add_argument("--min-dist", type=float, default=0.15)
    ap.add_argument("--n-jobs", type=int, default=16)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[INFO] read global QC", flush=True)
    qc = pd.read_csv(args.global_qc, sep="\t")
    qc = qc[qc["CB"].astype(str) != "CB"].copy()
    for c in ["reads_used", "covered_calls", "covered_bins", "modified_rate"]:
        qc[c] = pd.to_numeric(qc[c], errors="coerce")
    qc = qc.dropna().copy()

    print("[INFO] read region counts", flush=True)
    df = pd.read_csv(
        args.region_counts,
        sep="\t",
        compression="gzip",
        dtype={"CB": str, "region_id": str}
    )

    for c in ["modified_calls", "covered_calls", "modified_rate"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna().copy()

    qc_cells = set(qc["CB"])
    df = df[df["CB"].isin(qc_cells)].copy()

    print(f"[INFO] raw cell-region pairs={len(df)}", flush=True)

    cell_nreg = df.groupby("CB")["region_id"].nunique().rename("n_regions_observed")
    keep_cells = set(cell_nreg[cell_nreg >= args.min_cell_regions].index)

    df = df[df["CB"].isin(keep_cells)].copy()
    qc = qc[qc["CB"].isin(keep_cells)].copy()

    cells = sorted(df["CB"].unique())
    if args.max_cells > 0 and len(cells) > args.max_cells:
        rng = np.random.default_rng(args.seed)
        cells = sorted(rng.choice(cells, size=args.max_cells, replace=False))

    df = df[df["CB"].isin(set(cells))].copy()
    qc = qc[qc["CB"].isin(set(cells))].copy()
    qc = qc.set_index("CB").loc[cells].reset_index()

    print(f"[INFO] selected cells={len(cells)}", flush=True)

    region_stats = df.groupby("region_id")["modified_rate"].agg(["count", "mean", "var"]).reset_index()
    region_stats = region_stats[region_stats["count"] >= args.min_region_cells].copy()
    region_stats = region_stats.sort_values("var", ascending=False).head(args.top_variable_regions).copy()

    regions = list(region_stats["region_id"])
    df = df[df["region_id"].isin(set(regions))].copy()

    print(f"[INFO] selected regions={len(regions)}", flush=True)

    region_stats.to_csv(outdir / f"{args.prefix}.selected_regions.tsv", sep="\t", index=False)
    cell_nreg.loc[cells].to_csv(outdir / f"{args.prefix}.cell_region_coverage.tsv", sep="\t", header=True)

    print("[INFO] pivot matrix", flush=True)
    mat = df.pivot_table(
        index="CB",
        columns="region_id",
        values="modified_rate",
        aggfunc="mean"
    )
    mat = mat.reindex(index=cells, columns=regions)

    # region均值填补缺失，和文献里的imputation逻辑一致
    mat = mat.apply(lambda x: x.fillna(x.mean()), axis=0)
    X = mat.to_numpy(dtype=np.float32)

    if args.mode == "residual":
        g = qc["modified_rate"].to_numpy(dtype=float)
        X = logit(X) - logit(g[:, None])
        X = np.clip(X, -5, 5).astype(np.float32)
    elif args.mode == "logit_rate":
        X = logit(X).astype(np.float32)
    elif args.mode == "activity_5mc":
        X = np.log1p(1.0 / (X + 0.01)).astype(np.float32)

    print("[INFO] PCA", flush=True)
    Xs = StandardScaler(with_mean=True, with_std=True).fit_transform(X)
    pcs = PCA(n_components=args.n_pcs, random_state=args.seed).fit_transform(Xs)

    pcs_for_umap = pcs.copy()

    if args.regress_qc:
        print("[INFO] regress QC covariates from PCs", flush=True)
        qc_mat = np.column_stack([
            np.log10(qc["reads_used"].to_numpy(dtype=float) + 1.0),
            np.log10(qc["covered_calls"].to_numpy(dtype=float) + 1.0),
            np.log10(qc["covered_bins"].to_numpy(dtype=float) + 1.0),
            qc["modified_rate"].to_numpy(dtype=float)
        ])
        qc_mat = StandardScaler(with_mean=True, with_std=True).fit_transform(qc_mat)
        design = np.column_stack([np.ones(qc_mat.shape[0]), qc_mat])
        beta = np.linalg.lstsq(design, pcs_for_umap, rcond=None)[0]
        fitted_qc = qc_mat @ beta[1:, :]
        pcs_for_umap = pcs_for_umap - fitted_qc

    if args.drop_pc1:
        print("[INFO] drop PC1", flush=True)
        pcs_for_umap = pcs_for_umap[:, 1:]

    print("[INFO] UMAP", flush=True)
    emb = umap.UMAP(
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        metric="euclidean",
        random_state=None,
        n_jobs=args.n_jobs
    ).fit_transform(pcs_for_umap)

    print("[INFO] clustering", flush=True)
    cluster = KMeans(n_clusters=args.n_clusters, random_state=args.seed, n_init=20).fit_predict(pcs_for_umap)

    res = qc.copy()
    res["UMAP1"] = emb[:, 0]
    res["UMAP2"] = emb[:, 1]
    res["cluster"] = cluster

    res.to_csv(outdir / f"{args.prefix}.umap.tsv", sep="\t", index=False)

    corr = res[["UMAP1", "UMAP2", "reads_used", "covered_calls", "covered_bins", "modified_rate"]].corr(method="spearman")
    corr.to_csv(outdir / f"{args.prefix}.umap_qc_spearman_corr.tsv", sep="\t")

    cluster_counts = res["cluster"].value_counts().sort_index()
    cluster_counts.to_csv(outdir / f"{args.prefix}.cluster_counts.tsv", sep="\t", header=["n_cells"])

    def plot(col, discrete=False):
        plt.figure(figsize=(6, 5.5))
        if discrete:
            plt.scatter(
                res["UMAP1"], res["UMAP2"],
                c=res[col],
                s=1.5,
                alpha=0.70,
                cmap="tab10",
                linewidths=0
            )
        else:
            plt.scatter(
                res["UMAP1"], res["UMAP2"],
                c=res[col],
                s=1.5,
                alpha=0.70,
                linewidths=0
            )
        cb = plt.colorbar()
        cb.set_label(col)
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.title(f"{args.prefix} colored by {col}")
        plt.tight_layout()
        plt.savefig(outdir / f"{args.prefix}.umap_by_{col}.png", dpi=300)
        plt.close()

    plot("cluster", discrete=True)
    for col in ["modified_rate", "reads_used", "covered_calls", "covered_bins"]:
        plot(col)

    # density图适合论文展示，不会因为点太多糊成一片
    plt.figure(figsize=(6, 5.5))
    plt.hexbin(res["UMAP1"], res["UMAP2"], gridsize=80, mincnt=1)
    cb = plt.colorbar()
    cb.set_label("cell density")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(f"{args.prefix} cell density")
    plt.tight_layout()
    plt.savefig(outdir / f"{args.prefix}.umap_density_hexbin.png", dpi=300)
    plt.close()

    print("[INFO] done", flush=True)

if __name__ == "__main__":
    main()
