#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def pick_col(df, candidates):
    lower = {c.lower(): c for c in df.columns}
    for x in candidates:
        if x in df.columns:
            return x
        if x.lower() in lower:
            return lower[x.lower()]
    return None

def robust_z(x):
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if mad > 0:
        den = 1.4826 * mad
    else:
        den = np.nanstd(x)
        if den == 0 or np.isnan(den):
            den = 1.0
    return (x - med) / den

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", default="snTAPS-FL-1.strict_qc_core")
    ap.add_argument("--min-reads", type=float, default=100000)
    ap.add_argument("--qlo", type=float, default=0.20)
    ap.add_argument("--qhi", type=float, default=0.80)
    ap.add_argument("--targets", default="10000,5000,2000")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    qc = pd.read_csv(args.qc, sep="\t")
    print("[INFO] input QC:", args.qc)
    print("[INFO] input cells:", len(qc))
    print("[INFO] columns:", list(qc.columns))

    cb_col = pick_col(qc, ["CB", "cell", "barcode", "cell_barcode"])
    reads_col = pick_col(qc, ["reads_used", "reads", "read_count", "n_reads"])
    calls_col = pick_col(qc, ["covered_calls", "cytosine_calls", "covered_C_calls", "total_covered_calls"])
    bins_col = pick_col(qc, ["covered_bins", "bins_covered", "n_covered_bins"])
    rate_col = pick_col(qc, ["modified_rate", "methylation_rate", "overall_modified_rate"])

    if cb_col is None:
        raise ValueError("Cannot find cell barcode column")
    if reads_col is None:
        raise ValueError("Cannot find reads column")

    metric_cols = {
        "reads_used": reads_col,
        "covered_calls": calls_col,
        "covered_bins": bins_col,
        "modified_rate": rate_col,
    }
    metric_cols = {k: v for k, v in metric_cols.items() if v is not None}
    print("[INFO] metric columns:", metric_cols)

    mask = pd.Series(True, index=qc.index)
    mask &= qc[reads_col].astype(float) >= args.min_reads

    threshold_rows = []
    for name, col in metric_cols.items():
        vals = qc.loc[mask, col].astype(float)
        lo = vals.quantile(args.qlo)
        hi = vals.quantile(args.qhi)
        threshold_rows.append({
            "metric": name,
            "column": col,
            "qlo": args.qlo,
            "qhi": args.qhi,
            "lo_value": lo,
            "hi_value": hi,
        })
        mask &= qc[col].astype(float).between(lo, hi)

    sub = qc.loc[mask].copy()
    print("[INFO] cells after strict quantile filter:", len(sub))

    if len(sub) == 0:
        raise ValueError("No cells left after strict QC filter")

    z_list = []
    for name, col in metric_cols.items():
        x = sub[col].astype(float).to_numpy()
        if name != "modified_rate":
            x = np.log10(x + 1)
        z = robust_z(x)
        sub[f"qc_z_{name}"] = z
        z_list.append(z)

    zmat = np.vstack(z_list).T
    sub["qc_core_score"] = np.sqrt(np.nanmean(zmat ** 2, axis=1))
    sub = sub.sort_values("qc_core_score", ascending=True)

    all_file = outdir / f"{args.prefix}.all_after_qc.global_per_cb.tsv"
    sub.to_csv(all_file, sep="\t", index=False)

    pd.DataFrame(threshold_rows).to_csv(
        outdir / f"{args.prefix}.filter_thresholds.tsv",
        sep="\t",
        index=False
    )

    summary = []
    targets = [int(x) for x in args.targets.split(",") if x.strip()]
    for n in targets:
        take = sub.head(min(n, len(sub))).copy()
        tag = f"core{len(take)}"
        tsv = outdir / f"{args.prefix}.{tag}.global_per_cb.tsv"
        cells = outdir / f"{args.prefix}.{tag}.cells.txt"
        take.to_csv(tsv, sep="\t", index=False)
        take[[cb_col]].to_csv(cells, sep="\t", index=False, header=False)
        summary.append({
            "target": n,
            "output_cells": len(take),
            "global_qc": str(tsv),
            "cell_list": str(cells),
        })
        print("[INFO] wrote", tag, len(take), "cells")

    pd.DataFrame(summary).to_csv(
        outdir / f"{args.prefix}.core_cellsets.summary.tsv",
        sep="\t",
        index=False
    )

if __name__ == "__main__":
    main()
