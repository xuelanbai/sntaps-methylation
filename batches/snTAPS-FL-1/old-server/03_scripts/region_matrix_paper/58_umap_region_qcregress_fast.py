#!/usr/bin/env python3
import argparse
from pathlib import Path
import gzip
import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import umap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_table_auto(path, **kwargs):
    path = str(path)
    if path.endswith(".gz"):
        return pd.read_csv(path, sep="\t", compression="gzip", **kwargs)
    return pd.read_csv(path, sep="\t", **kwargs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region-counts", required=True)
    ap.add_argument("--cell-qc", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--mode", choices=["residual", "seurat_norm", "raw"], default="residual")
    ap.add_argument("--qc-keys", default="reads_used")
    ap.add_argument("--n-top-regions", type=int, default=3000)
    ap.add_argument("--min-region-cells", type=int, default=100)
    ap.add_argument("--min-regions-per-cell", type=int, default=20)
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--n-clusters", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--chunksize", type=int, default=2000000)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[INFO] read cell QC", flush=True)
    qc = read_table_auto(args.cell_qc)
    if "CB" not in qc.columns:
        raise ValueError("CB not found in cell QC")

    qc = qc.drop_duplicates("CB").copy()
    qc["CB"] = qc["CB"].astype(str)
    cb_list = qc["CB"].tolist()
    cb_set = set(cb_list)
    print(f"[INFO] cells from QC: {len(cb_list)}", flush=True)

    region_col = "region_id"
    cb_col = "CB"
    cov_col = "covered_calls"
    mod_col = "modified_calls"
    rate_col = "modified_rate"

    print("[INFO] first pass: region mean/variance among selected cells", flush=True)
    stats = {}

    usecols = None
    for chunk in pd.read_csv(args.region_counts, sep="\t", compression="gzip", chunksize=args.chunksize):
        if cb_col not in chunk.columns or region_col not in chunk.columns:
            raise ValueError(f"need columns {cb_col} and {region_col}; got {list(chunk.columns)}")
        sub = chunk[chunk[cb_col].astype(str).isin(cb_set)].copy()
        if sub.empty:
            continue

        sub[region_col] = sub[region_col].astype(str)
        sub[rate_col] = pd.to_numeric(sub[rate_col], errors="coerce")
        sub[cov_col] = pd.to_numeric(sub[cov_col], errors="coerce")
        sub = sub[np.isfinite(sub[rate_col]) & np.isfinite(sub[cov_col]) & (sub[cov_col] > 0)]
        if sub.empty:
            continue

        g = sub.groupby(region_col)[rate_col].agg(["count", "sum"])
        g2 = sub.assign(rate2=sub[rate_col] ** 2).groupby(region_col)["rate2"].sum()

        for rid in g.index:
            cnt = int(g.loc[rid, "count"])
            sm = float(g.loc[rid, "sum"])
            sm2 = float(g2.loc[rid])
            if rid not in stats:
                stats[rid] = [0, 0.0, 0.0]
            stats[rid][0] += cnt
            stats[rid][1] += sm
            stats[rid][2] += sm2

    rows = []
    for rid, (cnt, sm, sm2) in stats.items():
        if cnt >= args.min_region_cells:
            mean = sm / cnt
            var = max(sm2 / cnt - mean * mean, 0.0)
            rows.append((rid, cnt, mean, var))

    region_stats = pd.DataFrame(rows, columns=["region_id", "n_cells", "mean_rate", "var_rate"])
    if region_stats.empty:
        raise ValueError("no region left after min-region-cells filter")

    region_stats = region_stats.sort_values("var_rate", ascending=False)
    selected = region_stats.head(args.n_top_regions)["region_id"].astype(str).tolist()
    region_stats.to_csv(outdir / f"{args.prefix}.selected_region_stats.tsv", sep="\t", index=False)
    pd.Series(selected).to_csv(outdir / f"{args.prefix}.selected_regions.tsv", sep="\t", index=False, header=False)
    print(f"[INFO] selected regions: {len(selected)}", flush=True)

    cb_to_i = {cb:i for i, cb in enumerate(cb_list)}
    reg_to_j = {r:j for j, r in enumerate(selected)}

    X = np.full((len(cb_list), len(selected)), np.nan, dtype=np.float32)

    print("[INFO] second pass: build dense cell x selected-region matrix", flush=True)
    for chunk in pd.read_csv(args.region_counts, sep="\t", compression="gzip", chunksize=args.chunksize):
        sub = chunk[chunk[cb_col].astype(str).isin(cb_set)].copy()
        if sub.empty:
            continue
        sub[region_col] = sub[region_col].astype(str)
        sub = sub[sub[region_col].isin(reg_to_j)]
        if sub.empty:
            continue

        sub[rate_col] = pd.to_numeric(sub[rate_col], errors="coerce")
        sub[cov_col] = pd.to_numeric(sub[cov_col], errors="coerce")
        sub = sub[np.isfinite(sub[rate_col]) & np.isfinite(sub[cov_col]) & (sub[cov_col] > 0)]
        if sub.empty:
            continue

        if args.mode == "raw":
            sub["_value"] = sub[rate_col]
        elif args.mode == "residual":
            qtmp = qc.set_index("CB")["modified_rate"].to_dict() if "modified_rate" in qc.columns else {}
            sub["_cell_rate"] = sub[cb_col].astype(str).map(qtmp)
            sub["_value"] = sub[rate_col] - pd.to_numeric(sub["_cell_rate"], errors="coerce")
        elif args.mode == "seurat_norm":
            sub["_value"] = np.log1p(sub[rate_col] * 10000.0)

        sub = sub[np.isfinite(sub["_value"])]
        if sub.empty:
            continue

        piv = sub.pivot_table(index=cb_col, columns=region_col, values="_value", aggfunc="mean")
        row_idx = [cb_to_i[str(x)] for x in piv.index.astype(str)]
        col_idx = [reg_to_j[str(x)] for x in piv.columns.astype(str)]
        X[np.ix_(row_idx, col_idx)] = piv.to_numpy(dtype=np.float32)

    obs_regions = np.isfinite(X).sum(axis=1)
    keep_cell = obs_regions >= args.min_regions_per_cell
    print(f"[INFO] cells before selected-region coverage filter: {X.shape[0]}", flush=True)
    print(f"[INFO] cells after selected-region coverage filter: {keep_cell.sum()}", flush=True)

    if keep_cell.sum() < 50:
        raise ValueError("too few cells left after selected-region coverage filter")

    X = X[keep_cell, :]
    qc = qc.loc[keep_cell].reset_index(drop=True).copy()
    obs_regions = obs_regions[keep_cell]

    nan_fraction = float(np.isnan(X).mean())
    print(f"[INFO] matrix nan_fraction before fill: {nan_fraction:.4f}", flush=True)

    col_mean = np.nanmean(X, axis=0)
    bad = ~np.isfinite(col_mean)
    if bad.any():
        col_mean[bad] = 0
    inds = np.where(~np.isfinite(X))
    X[inds] = np.take(col_mean, inds[1])

    qc_keys = [x.strip() for x in args.qc_keys.split(",") if x.strip()]
    qc_keys = [x for x in qc_keys if x in qc.columns]
    print(f"[INFO] QC regress keys: {qc_keys}", flush=True)

    if len(qc_keys) > 0:
        Z = qc[qc_keys].copy()
        for c in qc_keys:
            Z[c] = pd.to_numeric(Z[c], errors="coerce")
            Z[c] = Z[c].fillna(Z[c].median())
            if Z[c].max() > 10:
                Z[c] = np.log10(Z[c] + 1)
        Z = StandardScaler(with_mean=True, with_std=True).fit_transform(Z.to_numpy(dtype=np.float64))
        design = np.column_stack([np.ones(Z.shape[0]), Z])
        beta = np.linalg.lstsq(design, X.astype(np.float64), rcond=None)[0]
        fitted = design[:, 1:] @ beta[1:, :]
        X = (X - fitted).astype(np.float32)

    print("[INFO] scale matrix after regression", flush=True)
    X = StandardScaler(with_mean=True, with_std=True).fit_transform(X).astype(np.float32)

    n_pcs = min(args.n_pcs, X.shape[0] - 1, X.shape[1] - 1)
    print(f"[INFO] PCA n_pcs={n_pcs}", flush=True)
    pcs = PCA(n_components=n_pcs, random_state=args.seed).fit_transform(X)

    print("[INFO] UMAP", flush=True)
    emb = umap.UMAP(
        n_neighbors=30,
        min_dist=0.3,
        metric="euclidean",
        random_state=args.seed,
        n_components=2,
        n_jobs=1
    ).fit_transform(pcs)

    print("[INFO] KMeans", flush=True)
    cluster = KMeans(n_clusters=args.n_clusters, random_state=args.seed, n_init=20).fit_predict(pcs)

    res = qc.copy()
    res["covered_selected_regions"] = obs_regions
    res["UMAP1"] = emb[:, 0]
    res["UMAP2"] = emb[:, 1]
    res["cluster"] = cluster
    res.to_csv(outdir / f"{args.prefix}.umap.tsv", sep="\t", index=False)

    res["cluster"].value_counts().sort_index().to_csv(
        outdir / f"{args.prefix}.cluster_counts.tsv",
        sep="\t",
        header=["n_cells"]
    )

    pd.DataFrame([{
        "n_cells": res.shape[0],
        "n_regions": X.shape[1],
        "nan_fraction_before_fill": nan_fraction,
        "mode": args.mode,
        "qc_keys": ",".join(qc_keys),
        "min_region_cells": args.min_region_cells,
        "min_regions_per_cell": args.min_regions_per_cell,
        "n_top_regions": args.n_top_regions
    }]).to_csv(outdir / f"{args.prefix}.matrix_qc.tsv", sep="\t", index=False)

    def scatter_color(col, suffix):
        plt.figure(figsize=(7, 6))
        if col == "cluster":
            sc = plt.scatter(res["UMAP1"], res["UMAP2"], c=res[col], s=4, cmap="tab10", linewidths=0)
            plt.colorbar(sc, label=col)
        else:
            val = pd.to_numeric(res[col], errors="coerce")
            sc = plt.scatter(res["UMAP1"], res["UMAP2"], c=val, s=4, cmap="viridis", linewidths=0)
            plt.colorbar(sc, label=col)
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.title(f"{args.prefix} colored by {col}")
        plt.tight_layout()
        plt.savefig(outdir / f"{args.prefix}.{suffix}.png", dpi=200)
        plt.close()

    scatter_color("cluster", "umap_by_cluster")
    for col in ["reads_used", "covered_calls", "covered_bins", "modified_rate", "covered_selected_regions"]:
        if col in res.columns:
            scatter_color(col, f"umap_by_{col}")

    print("[INFO] done", flush=True)


if __name__ == "__main__":
    main()
