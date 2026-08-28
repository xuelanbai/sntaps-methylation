#!/usr/bin/env python3
import re
import argparse
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_UMAP = "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"
DEFAULT_QC = "08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.global_per_cb.tsv"
DEFAULT_REGION = "08_qc_summary/final_single_sample/region_matrix_covfilter/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.cov10.region_counts.tsv.gz"
DEFAULT_MARKER = "08_qc_summary/final_single_sample/P7_marker_heatmap/P7.marker_promoter_activity_zscore_matrix.tsv"
DEFAULT_OUT = "08_qc_summary/final_single_sample/P7_cluster01_doublet_check"


def norm(x):
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())


def find_col(df, candidates, required=False, label="column"):
    cols = list(df.columns)

    for c in candidates:
        if c in cols:
            return c

    nmap = {norm(c): c for c in cols}
    for c in candidates:
        if norm(c) in nmap:
            return nmap[norm(c)]

    for c in candidates:
        k = norm(c)
        for col in cols:
            nk = norm(col)
            if k and (nk == k or nk.endswith(k) or k in nk):
                return col

    if required:
        raise ValueError(
            f"Cannot find {label}. Tried: {candidates}\n"
            f"Available columns:\n" + "\n".join(map(str, cols))
        )
    return None


def read_tsv(path, **kwargs):
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False, **kwargs)


def parse_targets(s):
    return [int(x) for x in re.split(r"[,; ]+", str(s).strip()) if x != ""]


def parse_cluster_id(x):
    s = str(x)
    m = re.search(r"(?:^|[\s_])C(?:luster)?\s*([0-9]+)\s*$", s, flags=re.I)
    if m:
        return int(m.group(1))
    m = re.search(r"([0-9]+)\s*$", s)
    if m:
        return int(m.group(1))
    return None


def clean_marker_name(x):
    s = str(x).strip()
    s = re.sub(r"\s+C(?:luster)?\s*[0-9]+\s*$", "", s, flags=re.I)
    s = re.sub(r"\s+cluster\s*[0-9]+\s*$", "", s, flags=re.I)
    return s.strip()


def marker_variants(x):
    base = clean_marker_name(x)
    vals = {base}
    if "|" in base:
        vals.add(base.split("|", 1)[0].strip())
    if " " in base:
        vals.add(base.split()[0].strip())
    return {v for v in vals if v and v.lower() not in {"nan", "none"}}


def load_marker_owner(marker_path):
    marker = pd.read_csv(marker_path, sep="\t", index_col=0)

    cluster_cols = {}
    for c in marker.columns:
        cid = parse_cluster_id(c)
        if cid is not None:
            cluster_cols[cid] = c

    tmp = defaultdict(set)
    records = []

    for raw_label, row in marker.iterrows():
        owner = None

        m = re.search(r"\s+C(?:luster)?\s*([0-9]+)\s*$", str(raw_label), flags=re.I)
        if m:
            owner = int(m.group(1))

        if owner is None and cluster_cols:
            vals = pd.to_numeric(row[list(cluster_cols.values())], errors="coerce")
            if vals.notna().any():
                best_col = vals.idxmax()
                for cid, col in cluster_cols.items():
                    if col == best_col:
                        owner = cid
                        break

        if owner is None:
            continue

        records.append({
            "raw_marker_label": raw_label,
            "marker_region": clean_marker_name(raw_label),
            "marker_cluster": owner,
        })

        for v in marker_variants(raw_label):
            tmp[v].add(owner)

    marker_owner = {k: list(v)[0] for k, v in tmp.items() if len(v) == 1}
    ambiguous = {k: sorted(v) for k, v in tmp.items() if len(v) > 1}

    return marker_owner, pd.DataFrame(records), ambiguous


def map_region_owner(region_series, marker_owner):
    s = region_series.astype(str)
    owner = s.map(marker_owner)

    miss = owner.isna()
    if miss.any():
        clean = s[miss].str.replace(r"\s+C(?:luster)?\s*[0-9]+\s*$", "", regex=True)
        owner.loc[miss] = clean.map(marker_owner)

    miss = owner.isna()
    if miss.any():
        prefix = s[miss].str.split("|", n=1, expand=True)[0].str.strip()
        owner.loc[miss] = prefix.map(marker_owner)

    miss = owner.isna()
    if miss.any():
        first = s[miss].str.split().str[0].str.strip()
        owner.loc[miss] = first.map(marker_owner)

    return owner


def detect_metric_cols(df):
    candidates = {
        "reads_used": [
            "reads_used", "unique_reads", "unique_read", "mapped_unique_reads",
            "dedup_reads", "deduplicated_reads", "valid_reads", "n_reads",
            "read_count", "reads"
        ],
        "covered_calls": [
            "covered_calls", "n_covered_calls", "total_calls",
            "valid_calls", "coverage_calls", "calls"
        ],
        "covered_bins": [
            "covered_bins", "n_covered_bins", "covered_regions",
            "n_regions", "regions", "bins"
        ],
        "modified_rate": [
            "modified_rate", "methylation_rate", "mc_rate",
            "mod_rate", "conversion_rate", "mean_modified_rate", "rate"
        ],
    }

    found = {}
    for metric, cands in candidates.items():
        col = find_col(df, cands, required=False, label=metric)
        if col is not None:
            found[metric] = col
    return found


def plot_umap_highlight(all_df, score_df, out_png):
    sus = set(score_df.loc[score_df["suspected_doublet_like"], "cell_barcode"].astype(str))
    targets = set(score_df["cell_barcode"].astype(str))

    p = all_df.copy()
    p["plot_group"] = "other"
    p.loc[p["cell_barcode"].astype(str).isin(targets), "plot_group"] = "cluster0_or_1"
    p.loc[p["cell_barcode"].astype(str).isin(sus), "plot_group"] = "suspected_doublet_like"

    fig, ax = plt.subplots(figsize=(7, 6))

    for group, color, size, alpha in [
        ("other", "#d9d9d9", 1, 0.12),
        ("cluster0_or_1", "#4c78a8", 2, 0.45),
        ("suspected_doublet_like", "#d62728", 5, 0.9),
    ]:
        sub = p[p["plot_group"] == group]
        if len(sub) == 0:
            continue
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=size, c=color, alpha=alpha, linewidths=0, label=group)

    ax.set_title("Cluster0/1 suspected doublet-like cells")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(frameon=False, markerscale=4)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umap", default=DEFAULT_UMAP)
    ap.add_argument("--qc", default=DEFAULT_QC)
    ap.add_argument("--region", default=DEFAULT_REGION)
    ap.add_argument("--marker", default=DEFAULT_MARKER)
    ap.add_argument("--outdir", default=DEFAULT_OUT)
    ap.add_argument("--targets", default="0,1")
    ap.add_argument("--chunksize", type=int, default=1000000)
    ap.add_argument("--top-z", type=float, default=1.0)
    ap.add_argument("--second-z", type=float, default=0.5)
    ap.add_argument("--margin-z", type=float, default=0.75)
    ap.add_argument("--qc-quantile", type=float, default=0.99)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("UMAP =", args.umap)
    print("QC =", args.qc)
    print("REGION_MATRIX =", args.region)
    print("MARKER =", args.marker)
    print("OUTDIR =", outdir)

    targets = parse_targets(args.targets)
    print("TARGET_CLUSTERS =", targets)

    umap = read_tsv(args.umap)
    if "Unnamed: 0" in umap.columns and find_col(umap, ["cell_barcode", "barcode", "CB", "cb", "cell"], required=False) is None:
        umap = umap.rename(columns={"Unnamed: 0": "cell_barcode"})

    cell_col = find_col(umap, ["cell_barcode", "barcode", "CB", "cb", "cell", "cell_id"], required=True, label="cell barcode in UMAP")
    cluster_col = find_col(umap, ["cluster", "leiden", "kmeans_cluster", "kmeans", "Cluster"], required=True, label="cluster")
    umap1_col = find_col(umap, ["UMAP1", "umap1", "UMAP_1", "x"], required=True, label="UMAP1")
    umap2_col = find_col(umap, ["UMAP2", "umap2", "UMAP_2", "y"], required=True, label="UMAP2")

    df = umap.rename(columns={
        cell_col: "cell_barcode",
        cluster_col: "orig_cluster",
        umap1_col: "UMAP1",
        umap2_col: "UMAP2",
    }).copy()

    df["cell_barcode"] = df["cell_barcode"].astype(str)
    df["orig_cluster"] = pd.to_numeric(df["orig_cluster"], errors="coerce").astype("Int64")

    # UMAP表里如果已有QC列，先读进来。
    for metric, col in detect_metric_cols(df).items():
        df[metric] = pd.to_numeric(df[col], errors="coerce")

    # 外部QC表：自动识别列名；没有reads_used也不退出。
    qc_path = Path(args.qc)
    if qc_path.exists():
        qc = read_tsv(qc_path)
        if "Unnamed: 0" in qc.columns and find_col(qc, ["cell_barcode", "barcode", "CB", "cb", "cell"], required=False) is None:
            qc = qc.rename(columns={"Unnamed: 0": "cell_barcode"})

        qc_cell = find_col(qc, ["cell_barcode", "barcode", "CB", "cb", "cell", "cell_id"], required=False, label="cell barcode in QC")
        if qc_cell is None:
            print("WARNING: cannot find cell barcode column in QC table; skip QC table.")
            print("QC columns:", list(qc.columns))
        else:
            qc = qc.rename(columns={qc_cell: "cell_barcode"})
            qc["cell_barcode"] = qc["cell_barcode"].astype(str)

            found = detect_metric_cols(qc)
            print("Detected QC metric columns:", found)

            keep = ["cell_barcode"] + list(found.values())
            qc2 = qc[keep].copy()
            qc2 = qc2.rename(columns={v: k for k, v in found.items()})

            for m in found.keys():
                qc2[m] = pd.to_numeric(qc2[m], errors="coerce")

            df = df.merge(qc2, on="cell_barcode", how="left", suffixes=("", "_qc"))

            for m in ["reads_used", "covered_calls", "covered_bins", "modified_rate"]:
                qcm = m + "_qc"
                if qcm in df.columns:
                    if m in df.columns:
                        df[m] = df[m].combine_first(df[qcm])
                    else:
                        df[m] = df[qcm]
                    df = df.drop(columns=[qcm])

    for m in ["reads_used", "covered_calls", "covered_bins", "modified_rate"]:
        if m not in df.columns:
            print(f"WARNING: {m} not found; skip this QC metric.")

    target_df = df[df["orig_cluster"].isin(targets)].copy()
    target_cells = set(target_df["cell_barcode"])
    print("Target cells =", len(target_cells))

    if len(target_cells) == 0:
        raise ValueError("No target cells found. Check --targets or cluster column.")

    # marker heatmap矩阵只提供marker列表；真正的per-cell marker activity从region matrix中重新计算。
    marker_owner, marker_records, ambiguous = load_marker_owner(args.marker)
    marker_records.to_csv(outdir / "marker_regions_used_for_signature.tsv", sep="\t", index=False)

    if ambiguous:
        pd.DataFrame(
            [{"marker_variant": k, "clusters": ",".join(map(str, v))} for k, v in ambiguous.items()]
        ).to_csv(outdir / "ambiguous_marker_variants_skipped.tsv", sep="\t", index=False)

    print("Usable marker variants =", len(marker_owner))

    header = read_tsv(args.region, nrows=5)

    reg_cell_col = find_col(header, ["cell_barcode", "barcode", "CB", "cb", "cell", "cell_id"], required=True, label="cell barcode in region matrix")
    region_col = find_col(header, ["region_id", "promoter", "gene_id", "gene", "region", "feature", "bin_id"], required=True, label="region id")

    rate_col = find_col(header, ["modified_rate", "methylation_rate", "mc_rate", "mod_rate", "rate"], required=False, label="modified_rate")
    mod_col = find_col(header, ["modified_calls", "n_modified", "modified", "methylated_calls", "mc"], required=False, label="modified calls")
    cov_col = find_col(header, ["covered_calls", "n_covered", "total_calls", "coverage", "cov", "calls"], required=False, label="covered calls")

    if rate_col is None and (mod_col is None or cov_col is None):
        raise ValueError(
            "Cannot compute marker activity. Need modified_rate, or modified_calls + covered_calls.\n"
            "Region matrix columns:\n" + "\n".join(map(str, header.columns))
        )

    print("Region matrix columns:")
    print("  cell =", reg_cell_col)
    print("  region =", region_col)
    print("  modified_rate =", rate_col)
    print("  modified_calls =", mod_col)
    print("  covered_calls =", cov_col)

    usecols = [reg_cell_col, region_col]
    if rate_col is not None:
        usecols.append(rate_col)
    else:
        usecols += [mod_col, cov_col]
    usecols = list(dict.fromkeys(usecols))

    sums = defaultdict(float)
    counts = defaultdict(float)

    total_rows = 0
    kept_rows = 0

    reader = pd.read_csv(
        args.region,
        sep="\t",
        compression="infer",
        usecols=usecols,
        chunksize=args.chunksize,
        low_memory=False,
    )

    for i, chunk in enumerate(reader, 1):
        total_rows += len(chunk)

        cell_s = chunk[reg_cell_col].astype(str)
        mask_cell = cell_s.isin(target_cells)

        if not mask_cell.any():
            continue

        sub = chunk.loc[mask_cell].copy()
        owner = map_region_owner(sub[region_col], marker_owner)
        mask_marker = owner.notna()

        if not mask_marker.any():
            continue

        sub = sub.loc[mask_marker].copy()
        owner = owner.loc[mask_marker].astype(int)
        sub["_marker_cluster"] = owner.values
        sub["_cell"] = sub[reg_cell_col].astype(str)

        if rate_col is not None:
            rate = pd.to_numeric(sub[rate_col], errors="coerce")
        else:
            mod = pd.to_numeric(sub[mod_col], errors="coerce")
            cov = pd.to_numeric(sub[cov_col], errors="coerce").replace(0, np.nan)
            rate = mod / cov

        # promoter hypomethylation activity：modified_rate越低，activity越高。
        sub["_activity"] = 1.0 - rate
        sub = sub.replace([np.inf, -np.inf], np.nan).dropna(subset=["_activity"])

        if len(sub) == 0:
            continue

        kept_rows += len(sub)

        g = sub.groupby(["_cell", "_marker_cluster"])["_activity"].agg(["sum", "count"])

        for (cell, clu), r in g.iterrows():
            sums[(cell, int(clu))] += float(r["sum"])
            counts[(cell, int(clu))] += float(r["count"])

        if i % 10 == 0:
            print(f"Processed chunks={i}; total_rows={total_rows}; marker_rows_kept={kept_rows}")

    marker_clusters = sorted(set(marker_owner.values()))

    sig = pd.DataFrame(index=sorted(target_cells))
    nobs = pd.DataFrame(index=sorted(target_cells))

    for c in marker_clusters:
        sig[f"marker_C{c}_activity"] = np.nan
        nobs[f"marker_C{c}_n"] = 0.0

    for (cell, c), s in sums.items():
        n = counts[(cell, c)]
        if n > 0:
            sig.loc[cell, f"marker_C{c}_activity"] = s / n
            nobs.loc[cell, f"marker_C{c}_n"] = n

    sig.index.name = "cell_barcode"
    nobs.index.name = "cell_barcode"

    sig.join(nobs).reset_index().to_csv(
        outdir / "cluster01_marker_activity_signature.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip"
    )

    sig_cols = [c for c in sig.columns if c.endswith("_activity")]

    Z = sig.copy()
    for c in sig_cols:
        mu = np.nanmean(Z[c].values)
        sd = np.nanstd(Z[c].values)
        if not np.isfinite(sd) or sd == 0:
            Z[c] = np.nan
        else:
            Z[c] = (Z[c] - mu) / sd

    arr = Z[sig_cols].to_numpy(dtype=float)
    arr2 = np.where(np.isfinite(arr), arr, -np.inf)
    order = np.argsort(arr2, axis=1)[:, ::-1]

    top1_idx = order[:, 0]
    top2_idx = order[:, 1] if len(sig_cols) > 1 else order[:, 0]

    top1_z = arr2[np.arange(arr2.shape[0]), top1_idx]
    top2_z = arr2[np.arange(arr2.shape[0]), top2_idx]

    def sig_cluster(col):
        m = re.search(r"marker_C([0-9]+)_activity", col)
        return int(m.group(1)) if m else -1

    top_df = pd.DataFrame(index=sig.index)
    top_df["marker_groups_observed"] = sig.notna().sum(axis=1)
    top_df["top1_marker_cluster"] = [sig_cluster(sig_cols[i]) for i in top1_idx]
    top_df["top2_marker_cluster"] = [sig_cluster(sig_cols[i]) for i in top2_idx]
    top_df["top1_marker_z"] = top1_z
    top_df["top2_marker_z"] = top2_z
    top_df["marker_margin_z"] = top_df["top1_marker_z"] - top_df["top2_marker_z"]

    scores = (
        target_df.set_index("cell_barcode")
        .join(sig, how="left")
        .join(nobs, how="left")
        .join(top_df, how="left")
        .reset_index()
    )

    scores["marker_mixed_like"] = (
        (scores["marker_groups_observed"] >= 2) &
        (scores["top1_marker_z"] >= args.top_z) &
        (scores["top2_marker_z"] >= args.second_z) &
        (scores["marker_margin_z"] <= args.margin_z)
    )

    qc_metrics = [m for m in ["reads_used", "covered_calls", "covered_bins"] if m in df.columns and df[m].notna().any()]

    for m in qc_metrics:
        cutoff = df[m].quantile(args.qc_quantile)
        scores[f"{m}_q{args.qc_quantile}_cutoff"] = cutoff
        scores[f"{m}_high"] = scores[m] >= cutoff

    high_cols = [f"{m}_high" for m in qc_metrics]

    if high_cols:
        scores["qc_high_metric_n"] = scores[high_cols].sum(axis=1)
        scores["qc_high_like"] = scores["qc_high_metric_n"] >= 2
    else:
        scores["qc_high_metric_n"] = 0
        scores["qc_high_like"] = False

    scores["suspected_doublet_like"] = scores["marker_mixed_like"] | scores["qc_high_like"]

    scores.to_csv(
        outdir / "cluster01_doublet_marker_qc_scores.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip"
    )

    suspected = scores[scores["suspected_doublet_like"]].copy()
    suspected.to_csv(
        outdir / "suspected_doublet_like_cluster01_cells.tsv",
        sep="\t",
        index=False
    )

    summary = []
    for cl in targets:
        sub = scores[scores["orig_cluster"].astype("Int64") == cl]
        summary.append({
            "cluster": cl,
            "n_cells": len(sub),
            "marker_mixed_like_n": int(sub["marker_mixed_like"].sum()),
            "qc_high_like_n": int(sub["qc_high_like"].sum()),
            "suspected_doublet_like_n": int(sub["suspected_doublet_like"].sum()),
            "suspected_doublet_like_frac": float(sub["suspected_doublet_like"].mean()) if len(sub) else np.nan,
            "median_reads_used": float(sub["reads_used"].median()) if "reads_used" in sub.columns else np.nan,
            "median_covered_calls": float(sub["covered_calls"].median()) if "covered_calls" in sub.columns else np.nan,
            "median_covered_bins": float(sub["covered_bins"].median()) if "covered_bins" in sub.columns else np.nan,
            "median_modified_rate": float(sub["modified_rate"].median()) if "modified_rate" in sub.columns else np.nan,
        })

    summary = pd.DataFrame(summary)
    summary.to_csv(outdir / "cluster01_doublet_check_summary.tsv", sep="\t", index=False)

    if len(suspected) > 0:
        pair = (
            suspected.groupby(["orig_cluster", "top1_marker_cluster", "top2_marker_cluster"])
            .size()
            .reset_index(name="n_cells")
            .sort_values(["orig_cluster", "n_cells"], ascending=[True, False])
        )
        pair.to_csv(outdir / "suspected_doublet_marker_pair_summary.tsv", sep="\t", index=False)

    plot_umap_highlight(df, scores, outdir / "cluster01_suspected_doublet_like_umap.png")

    print("\nDone.")
    print(summary.to_string(index=False))
    print("\nOutputs:")
    for p in sorted(outdir.glob("*")):
        print(p)


if __name__ == "__main__":
    main()
