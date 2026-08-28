#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
import pandas as pd
import numpy as np

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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--min-reads", type=float, default=100000)
    ap.add_argument("--filter-mode", choices=["reads_only", "auto"], default="reads_only")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    qc = pd.read_csv(args.qc, sep="\t")
    qc.columns = [str(c).strip() for c in qc.columns]
    cols = list(qc.columns)

    cb_col = find_col(
        cols,
        ["CB", "cell", "barcode", "cell_barcode", "cell_id"],
        r"^(CB|cell|barcode|cell_barcode|cell_id)$"
    )
    if cb_col is None:
        raise ValueError("Cannot find cell barcode column. Columns are:\n" + "\n".join(cols))

    reads_col = find_col(
        cols,
        ["reads_used", "reads", "read_count", "n_reads", "total_reads",
         "reads_per_cell", "mapped_reads", "unique_reads", "fragments", "n_fragments"],
        r"read|fragment"
    )
    if reads_col is None:
        raise ValueError("Cannot find reads column. Columns are:\n" + "\n".join(cols))

    metric_cols = [reads_col]

    if args.filter_mode == "auto":
        cyt_col = find_col(
            cols,
            ["covered_calls", "cytosine_calls", "cytosine_number", "n_cytosine",
             "total_covered_calls", "covered_cytosines"],
            r"covered.*call|cytosine"
        )
        bin_col = find_col(
            cols,
            ["covered_bins", "bins", "n_bins", "covered_100kb_bins", "bin_number"],
            r"covered.*bin|100.*bin|bin_number"
        )
        for c in [cyt_col, bin_col]:
            if c is not None and c not in metric_cols:
                metric_cols.append(c)

    print("[INFO] CB column:", cb_col, flush=True)
    print("[INFO] reads column:", reads_col, flush=True)
    print("[INFO] filter columns:", metric_cols, flush=True)

    qc[cb_col] = qc[cb_col].astype(str)
    mask = pd.Series(True, index=qc.index)

    summary = []
    for col in metric_cols:
        x = pd.to_numeric(qc[col], errors="coerce")
        q10 = float(np.nanquantile(x, 0.10))
        q50 = float(np.nanquantile(x, 0.50))
        q95 = float(np.nanquantile(x, 0.95))

        if col == reads_col:
            lower = max(q10, args.min_reads)
        else:
            lower = q10

        before = int(mask.sum())
        cur = x.notna() & (x >= lower) & (x <= q95)
        mask = mask & cur
        after = int(mask.sum())

        summary.append({
            "metric": col,
            "q10": q10,
            "median": q50,
            "q95": q95,
            "lower_used": lower,
            "upper_used": q95,
            "n_before_step": before,
            "n_after_step": after
        })

    qc["keep_strict"] = mask.values
    kept = qc.loc[mask].copy()

    qc.to_csv(outdir / f"{args.prefix}.all_cells_with_filter.tsv", sep="\t", index=False)
    kept.to_csv(outdir / f"{args.prefix}.strict_qc_cells.tsv", sep="\t", index=False)
    kept[[cb_col]].to_csv(outdir / f"{args.prefix}.strict_cells.txt", sep="\t", index=False, header=False)

    pd.DataFrame(summary).to_csv(outdir / f"{args.prefix}.filter_thresholds.tsv", sep="\t", index=False)

    print("[INFO] total cells:", len(qc), flush=True)
    print("[INFO] kept cells:", len(kept), flush=True)
    print("[INFO] output:", outdir, flush=True)

if __name__ == "__main__":
    main()
