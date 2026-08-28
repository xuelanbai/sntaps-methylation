#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
import numpy as np
import pandas as pd

def find_col(cols, exact_list=None, regex=None):
    exact_list = exact_list or []
    lower_map = {c.lower(): c for c in cols}
    for x in exact_list:
        if x.lower() in lower_map:
            return lower_map[x.lower()]
    if regex:
        for c in cols:
            if re.search(regex, c, flags=re.I):
                return c
    return None

def region_key_from_chunk(chunk, region_col, chrom_col, start_col, end_col):
    if region_col is not None:
        return chunk[region_col].astype(str)
    return (
        chunk[chrom_col].astype(str) + ":" +
        chunk[start_col].astype(str) + "-" +
        chunk[end_col].astype(str)
    )

def get_value(chunk, mod_col, cov_col, rate_col):
    if mod_col is not None and cov_col is not None:
        mod = pd.to_numeric(chunk[mod_col], errors="coerce")
        cov = pd.to_numeric(chunk[cov_col], errors="coerce")
        v = mod / cov.replace(0, np.nan)
        return v.astype("float32")
    else:
        return pd.to_numeric(chunk[rate_col], errors="coerce").astype("float32")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region-counts", required=True)
    ap.add_argument("--cell-qc", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)

    ap.add_argument("--n-top-regions", type=int, default=3000)
    ap.add_argument("--min-region-cells", type=int, default=100)
    ap.add_argument("--min-regions-per-cell", type=int, default=100)

    ap.add_argument("--n-pcs", type=int, default=50)
    ap.add_argument("--drop-pc1", action="store_true")
    ap.add_argument("--n-components", type=int, default=2)
    ap.add_argument("--n-neighbors", type=int, default=30)
    ap.add_argument("--min-dist", type=float, default=0.30)
    ap.add_argument("--n-clusters", type=int, default=8)
    ap.add_argument("--seed", type=int, default=-1)
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--chunksize", type=int, default=1000000)
    ap.add_argument("--point-size", type=float, default=2.0)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[INFO] read cell QC", flush=True)
    qc = pd.read_csv(args.cell_qc, sep="\t")
    qc.columns = [str(c).strip() for c in qc.columns]
    qc_cols = list(qc.columns)

    qc_cb_col = find_col(qc_cols, ["CB", "cell", "barcode", "cell_barcode", "cell_id"], r"^(CB|cell|barcode|cell_barcode|cell_id)$")
    if qc_cb_col is None:
        raise ValueError("Cannot find CB column in cell QC")

    qc[qc_cb_col] = qc[qc_cb_col].astype(str)
    qc = qc.drop_duplicates(qc_cb_col).reset_index(drop=True)
    cell_order = qc[qc_cb_col].tolist()
    cell_set = set(cell_order)
    print("[INFO] filtered cells input:", len(cell_order), flush=True)

    print("[INFO] inspect region_counts header", flush=True)
    header = pd.read_csv(args.region_counts, sep="\t", nrows=5)
    header.columns = [str(c).strip() for c in header.columns]
    cols = list(header.columns)
    print("[INFO] region_counts columns:", cols, flush=True)

    cb_col = find_col(cols, ["CB", "cell", "barcode", "cell_barcode", "cell_id"], r"^(CB|cell|barcode|cell_barcode|cell_id)$")
    if cb_col is None:
        raise ValueError("Cannot find CB column in region_counts")

    region_col = find_col(cols, ["region", "region_id", "promoter_id", "gene_id", "bin_id", "feature_id"], r"region|promoter|gene_id|bin_id|feature_id")
    chrom_col = find_col(cols, ["chrom", "chr", "seqnames"])
    start_col = find_col(cols, ["start", "region_start"])
    end_col = find_col(cols, ["end", "region_end"])

    if region_col is None:
        if chrom_col is None or start_col is None or end_col is None:
            raise ValueError("Cannot find region column, nor chrom/start/end columns")
        print("[INFO] use chrom:start-end as region ID", flush=True)
    else:
        print("[INFO] region column:", region_col, flush=True)

    cov_col = find_col(cols, ["covered_calls", "coverage", "cov", "total_calls", "n_calls", "cytosine_calls"], r"covered.*call|coverage|cytosine")
    mod_col = find_col(cols, ["modified_calls", "methylated_calls", "mod_calls", "meth_calls"], r"modified.*call|methylated.*call|meth.*call")
    rate_col = find_col(cols, ["modified_rate", "methylation_rate", "meth_rate", "ratio", "value"], r"modified_rate|methylation_rate|meth_rate|ratio|value")

    if not ((mod_col is not None and cov_col is not None) or rate_col is not None):
        raise ValueError("Cannot infer methylation value columns")

    print("[INFO] CB column:", cb_col, flush=True)
    print("[INFO] cov column:", cov_col, flush=True)
    print("[INFO] mod column:", mod_col, flush=True)
    print("[INFO] rate column:", rate_col, flush=True)

    usecols = [cb_col]
    if region_col is not None:
        usecols.append(region_col)
    else:
        usecols += [chrom_col, start_col, end_col]
    if mod_col is not None:
        usecols.append(mod_col)
    if cov_col is not None:
        usecols.append(cov_col)
    if rate_col is not None:
        usecols.append(rate_col)
    usecols = list(dict.fromkeys(usecols))

    print("[INFO] first pass: region variance", flush=True)
    stat = None
    total_rows = 0
    kept_rows = 0

    for chunk in pd.read_csv(args.region_counts, sep="\t", usecols=usecols, chunksize=args.chunksize):
        total_rows += len(chunk)
        chunk[cb_col] = chunk[cb_col].astype(str)
        chunk = chunk[chunk[cb_col].isin(cell_set)]
        kept_rows += len(chunk)
        if len(chunk) == 0:
            continue

        reg = region_key_from_chunk(chunk, region_col, chrom_col, start_col, end_col)
        v = get_value(chunk, mod_col, cov_col, rate_col)

        tmp = pd.DataFrame({"region": reg.values, "v": v.values})
        tmp = tmp.replace([np.inf, -np.inf], np.nan).dropna()
        if len(tmp) == 0:
            continue
        tmp["v2"] = tmp["v"] * tmp["v"]

        g = tmp.groupby("region", sort=False).agg(
            n=("v", "count"),
            sum=("v", "sum"),
            sumsq=("v2", "sum")
        )

        if stat is None:
            stat = g
        else:
            stat = stat.add(g, fill_value=0)

    print("[INFO] total rows:", total_rows, flush=True)
    print("[INFO] rows from filtered cells:", kept_rows, flush=True)

    if stat is None or len(stat) == 0:
        raise ValueError("No region-count rows found for filtered cells")

    stat["mean"] = stat["sum"] / stat["n"]
    stat["var"] = stat["sumsq"] / stat["n"] - stat["mean"] * stat["mean"]
    stat["var"] = stat["var"].clip(lower=0)

    eligible = stat[stat["n"] >= args.min_region_cells].copy()
    if len(eligible) < 500:
        relaxed = max(10, args.min_region_cells // 5)
        print(f"[WARN] too few eligible regions; relax min_region_cells to {relaxed}", flush=True)
        eligible = stat[stat["n"] >= relaxed].copy()

    selected = eligible.sort_values("var", ascending=False).head(args.n_top_regions).copy()
    selected = selected.reset_index().rename(columns={"index": "region"})
    selected.to_csv(outdir / f"{args.prefix}.selected_regions.tsv", sep="\t", index=False)

    selected_regions = selected["region"].astype(str).tolist()
    region_to_j = {r: i for i, r in enumerate(selected_regions)}

    print("[INFO] selected regions:", len(selected_regions), flush=True)

    print("[INFO] second pass: build dense selected matrix", flush=True)
    cell_to_i = {c: i for i, c in enumerate(cell_order)}
    X = np.full((len(cell_order), len(selected_regions)), np.nan, dtype=np.float32)

    selected_set = set(selected_regions)

    for chunk in pd.read_csv(args.region_counts, sep="\t", usecols=usecols, chunksize=args.chunksize):
        chunk[cb_col] = chunk[cb_col].astype(str)
        chunk = chunk[chunk[cb_col].isin(cell_set)]
        if len(chunk) == 0:
            continue

        reg = region_key_from_chunk(chunk, region_col, chrom_col, start_col, end_col)
        keep = reg.isin(selected_set)
        if keep.sum() == 0:
            continue

        chunk = chunk.loc[keep].copy()
        reg = reg.loc[keep]
        v = get_value(chunk, mod_col, cov_col, rate_col)

        ii = chunk[cb_col].map(cell_to_i)
        jj = reg.map(region_to_j)

        valid = ii.notna() & jj.notna() & v.notna() & np.isfinite(v)
        if valid.sum() == 0:
            continue

        X[ii[valid].astype(int).values, jj[valid].astype(int).values] = v[valid].values.astype(np.float32)

    obs_per_cell = np.sum(~np.isnan(X), axis=1)
    cell_cov = pd.DataFrame({
        "CB": cell_order,
        "covered_selected_regions": obs_per_cell
    })
    cell_cov.to_csv(outdir / f"{args.prefix}.cell_region_coverage.tsv", sep="\t", index=False)

    keep_cell = obs_per_cell >= args.min_regions_per_cell
    print("[INFO] cells before selected-region coverage filter:", len(cell_order), flush=True)
    print("[INFO] cells after selected-region coverage filter:", int(keep_cell.sum()), flush=True)

    X = X[keep_cell, :]
    qc = qc.loc[keep_cell].reset_index(drop=True)
    cell_order = [c for c, k in zip(cell_order, keep_cell) if k]
    obs_per_cell = obs_per_cell[keep_cell]

    if X.shape[0] < 50:
        raise ValueError("Too few cells left after selected-region coverage filter")

    print("[INFO] impute missing by region mean", flush=True)
    region_mean = np.nanmean(X, axis=0)
    global_mean = np.nanmean(region_mean)
    region_mean = np.where(np.isfinite(region_mean), region_mean, global_mean)

    inds = np.where(np.isnan(X))
    X[inds] = np.take(region_mean, inds[1])

    print("[INFO] scale matrix", flush=True)
    X = X.astype(np.float32)
    X -= X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, keepdims=True)
    sd[sd == 0] = 1
    X /= sd

    print("[INFO] PCA", flush=True)
    from sklearn.decomposition import PCA
    n_pcs = min(args.n_pcs, X.shape[0] - 1, X.shape[1] - 1)
    pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=None if args.seed < 0 else args.seed)
    pcs = pca.fit_transform(X)

    pd.DataFrame({
        "PC": np.arange(1, len(pca.explained_variance_ratio_) + 1),
        "explained_variance_ratio": pca.explained_variance_ratio_
    }).to_csv(outdir / f"{args.prefix}.pca_explained_variance.tsv", sep="\t", index=False)

    pcs_for_umap = pcs.copy()

    if args.drop_pc1 and pcs_for_umap.shape[1] > 2:
        print("[INFO] drop PC1", flush=True)
        pcs_for_umap = pcs_for_umap[:, 1:]

    print("[INFO] UMAP", flush=True)
    import umap
    reducer = umap.UMAP(
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        metric="euclidean",
        n_components=args.n_components,
        random_state=None if args.seed < 0 else args.seed,
        n_jobs=args.n_jobs
    )
    emb = reducer.fit_transform(pcs_for_umap)

    print("[INFO] clustering", flush=True)
    cluster = None
    try:
        import scanpy as sc
        import anndata as ad
        adata = ad.AnnData(pcs_for_umap)
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep="X")
        sc.tl.leiden(adata, resolution=0.6, random_state=0 if args.seed < 0 else args.seed)
        cluster = adata.obs["leiden"].astype(str).values
        cluster_method = "leiden"
    except Exception as e:
        print("[WARN] Leiden failed; fallback to MiniBatchKMeans:", repr(e), flush=True)
        from sklearn.cluster import MiniBatchKMeans
        cluster = MiniBatchKMeans(
            n_clusters=args.n_clusters,
            random_state=None if args.seed < 0 else args.seed,
            n_init=20,
            batch_size=4096
        ).fit_predict(pcs_for_umap).astype(str)
        cluster_method = "MiniBatchKMeans"

    res = qc.copy()
    res["CB"] = cell_order
    for i in range(args.n_components):
        res[f"UMAP{i+1}"] = emb[:, i]
    res["cluster"] = cluster
    res["cluster_method"] = cluster_method
    res["covered_selected_regions"] = obs_per_cell

    res.to_csv(outdir / f"{args.prefix}.umap.tsv", sep="\t", index=False)

    counts = res["cluster"].value_counts().sort_index()
    counts.to_csv(outdir / f"{args.prefix}.cluster_counts.tsv", sep="\t", header=["n_cells"])

    corr_cols = [c for c in res.columns if c.startswith("UMAP")]
    for c in res.columns:
        if re.search(r"read|covered|cytosine|bin|modified|meth", c, flags=re.I):
            if pd.api.types.is_numeric_dtype(res[c]):
                corr_cols.append(c)
    corr_cols = list(dict.fromkeys(corr_cols))
    if len(corr_cols) > 2:
        res[corr_cols].corr(method="spearman").to_csv(
            outdir / f"{args.prefix}.umap_qc_spearman_corr.tsv",
            sep="\t"
        )

    print("[INFO] plotting", flush=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    codes = pd.Categorical(res["cluster"]).codes

    plt.figure(figsize=(7, 6))
    plt.scatter(res["UMAP1"], res["UMAP2"], c=codes, s=args.point_size, alpha=0.75, linewidths=0, rasterized=True)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(f"{args.prefix} cluster")
    plt.tight_layout()
    plt.savefig(outdir / f"{args.prefix}.umap_by_cluster.png", dpi=300)
    plt.close()

    plt.figure(figsize=(7, 6))
    plt.hexbin(res["UMAP1"], res["UMAP2"], gridsize=80, mincnt=1)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(f"{args.prefix} density")
    plt.tight_layout()
    plt.savefig(outdir / f"{args.prefix}.umap_density_hexbin.png", dpi=300)
    plt.close()

    for color_col in ["covered_selected_regions"]:
        plt.figure(figsize=(7, 6))
        plt.scatter(res["UMAP1"], res["UMAP2"], c=res[color_col], s=args.point_size, alpha=0.75, linewidths=0, rasterized=True)
        plt.colorbar(label=color_col)
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.title(f"{args.prefix} by {color_col}")
        plt.tight_layout()
        plt.savefig(outdir / f"{args.prefix}.umap_by_{color_col}.png", dpi=300)
        plt.close()

    if args.n_components >= 3:
        try:
            import plotly.express as px
            html_dir = outdir / "plotly_html"
            html_dir.mkdir(exist_ok=True)
            fig = px.scatter_3d(
                res,
                x="UMAP1",
                y="UMAP2",
                z="UMAP3",
                color="cluster",
                hover_data=["CB", "cluster", "covered_selected_regions"],
                title=f"{args.prefix} 3D UMAP"
            )
            fig.update_traces(marker=dict(size=2, opacity=0.75))
            fig.write_html(html_dir / f"{args.prefix}.3d_umap_by_cluster.html")
        except Exception as e:
            print("[WARN] plotly 3D html failed:", repr(e), flush=True)

    print("[INFO] done:", outdir, flush=True)

if __name__ == "__main__":
    main()
