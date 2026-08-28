#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path
from collections import defaultdict

import pandas as pd


def pick_col(cols, candidates, label):
    for c in candidates:
        if c in cols:
            return c
    raise ValueError(f"Cannot find {label} column. Header={cols}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--umap", required=True, help="UMAP tsv with CB and cluster columns")
    parser.add_argument("--site", required=True, help="Merged CpG site conversion tsv.gz")
    parser.add_argument("--chrom", required=True, help="Chromosome to process, e.g. 1 or chr1")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--mincov", type=int, default=5)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. read cluster assignment
    umap = pd.read_csv(args.umap, sep="\t")
    if "CB" not in umap.columns or "cluster" not in umap.columns:
        raise ValueError(f"UMAP file must contain CB and cluster columns. Found: {umap.columns.tolist()}")

    umap = umap[["CB", "cluster"]].dropna().copy()
    umap["cluster"] = umap["cluster"].astype(int)

    cb_to_cluster = dict(zip(umap["CB"], umap["cluster"]))
    clusters = sorted(umap["cluster"].unique().tolist())
    selected_cbs = set(cb_to_cluster.keys())

    print(f"[INFO] selected cells: {len(selected_cbs)}", flush=True)
    print(f"[INFO] clusters: {clusters}", flush=True)
    print(f"[INFO] target chrom: {args.chrom}", flush=True)

    # 2. detect site file columns
    with gzip.open(args.site, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")

    chrom_col = pick_col(header, ["chrom", "chr", "#chrom", "chromosome"], "chrom")
    pos_col   = pick_col(header, ["pos", "position", "start"], "position")
    cb_col    = pick_col(header, ["CB", "cb", "barcode", "cell_barcode", "cell"], "CB")
    mod_col   = pick_col(header, ["modified_calls", "modified_count", "meth_calls", "methylated_calls", "n_mod"], "modified_calls")
    cov_col   = pick_col(header, ["covered_calls", "coverage", "cov", "total_calls", "n_cov"], "covered_calls")

    usecols = [chrom_col, pos_col, cb_col, mod_col, cov_col]
    print(f"[INFO] site columns used: {usecols}", flush=True)

    # 3. aggregate by cluster + position
    agg = {c: defaultdict(lambda: [0, 0]) for c in clusters}

    chunk_id = 0
    reader = pd.read_csv(args.site, sep="\t", usecols=usecols, chunksize=1_000_000)
    for chunk in reader:
        chunk_id += 1
        if chunk_id % 10 == 0:
            print(f"[INFO] reading chunk {chunk_id}", flush=True)

        chunk[chrom_col] = chunk[chrom_col].astype(str)
        sub = chunk[chunk[chrom_col] == str(args.chrom)].copy()
        if sub.empty:
            continue

        sub = sub[sub[cb_col].isin(selected_cbs)].copy()
        if sub.empty:
            continue

        sub["cluster"] = sub[cb_col].map(cb_to_cluster)
        sub = sub.dropna(subset=["cluster"])
        if sub.empty:
            continue

        sub["cluster"] = sub["cluster"].astype(int)
        sub[pos_col] = sub[pos_col].astype(int)
        sub[mod_col] = sub[mod_col].astype(float)
        sub[cov_col] = sub[cov_col].astype(float)

        grp = (
            sub.groupby(["cluster", pos_col], as_index=False)[[mod_col, cov_col]]
            .sum()
        )

        for row in grp.itertuples(index=False):
            cluster = int(getattr(row, "cluster"))
            pos = int(getattr(row, pos_col))
            mod = float(getattr(row, mod_col))
            cov = float(getattr(row, cov_col))
            agg[cluster][pos][0] += mod
            agg[cluster][pos][1] += cov

    # 4. write bedGraphs
    for c in clusters:
        meth_bg = outdir / f"C{c}.{args.chrom}.methyl.bedGraph"
        hypo_bg = outdir / f"C{c}.{args.chrom}.hypo.bedGraph"
        cov_bg  = outdir / f"C{c}.{args.chrom}.coverage.bedGraph"

        n_written = 0
        with open(meth_bg, "w") as fm, open(hypo_bg, "w") as fh, open(cov_bg, "w") as fc:
            for pos in sorted(agg[c].keys()):
                mod, cov = agg[c][pos]
                if cov < args.mincov:
                    continue
                start = pos - 1
                end = pos
                meth = mod / cov
                hypo = 1.0 - meth

                fm.write(f"{args.chrom}\t{start}\t{end}\t{meth:.6f}\n")
                fh.write(f"{args.chrom}\t{start}\t{end}\t{hypo:.6f}\n")
                fc.write(f"{args.chrom}\t{start}\t{end}\t{cov:.2f}\n")
                n_written += 1

        print(f"[INFO] C{c} wrote {n_written} sites on {args.chrom}", flush=True)


if __name__ == "__main__":
    main()