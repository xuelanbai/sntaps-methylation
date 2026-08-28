#!/usr/bin/env python3
import argparse
import gzip
import os
import subprocess
from pathlib import Path
import pandas as pd
import numpy as np

def find_col(cols, names):
    lower = {c.lower(): c for c in cols}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    for c in cols:
        cl = c.lower()
        for n in names:
            if n.lower() in cl:
                return c
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--genome-size", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-cov", type=int, default=5)
    ap.add_argument("--chunk-size", type=int, default=10000000)
    ap.add_argument("--bedGraphToBigWig", default="bedGraphToBigWig")
    args = ap.parse_args()

    infile = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    prefix = infile.name.replace(".tsv.gz", "").replace(".tsv", "")
    print(f"[INFO] input = {infile}", flush=True)

    # 先读表头
    with gzip.open(infile, "rt") if str(infile).endswith(".gz") else open(infile) as f:
        header = f.readline().rstrip("\n").split("\t")

    chrom_col = find_col(header, ["chrom", "chr", "chromosome", "seqnames", "contig"])
    start_col = find_col(header, ["start", "chromStart"])
    end_col = find_col(header, ["end", "chromEnd"])
    pos_col = find_col(header, ["pos", "position", "site"])

    mod_col = find_col(header, ["modified_calls", "methylated_calls", "meth_calls", "mod_calls", "modified"])
    cov_col = find_col(header, ["covered_calls", "coverage", "cov", "depth", "total_calls", "total_covered_calls"])
    rate_col = find_col(header, ["modified_rate", "methylation_rate", "meth_rate", "ratio", "conversion_rate"])

    print("[INFO] detected columns:", flush=True)
    print("  chrom =", chrom_col, flush=True)
    print("  start =", start_col, flush=True)
    print("  end   =", end_col, flush=True)
    print("  pos   =", pos_col, flush=True)
    print("  mod   =", mod_col, flush=True)
    print("  cov   =", cov_col, flush=True)
    print("  rate  =", rate_col, flush=True)

    if chrom_col is None:
        raise ValueError("Cannot detect chromosome column")
    if cov_col is None:
        raise ValueError("Cannot detect coverage column")
    if rate_col is None and (mod_col is None or cov_col is None):
        raise ValueError("Need either rate column, or modified_calls + covered_calls columns")
    if pos_col is None and start_col is None:
        raise ValueError("Need pos column or start column")

    usecols = [chrom_col, cov_col]
    if pos_col is not None:
        usecols.append(pos_col)
    else:
        usecols.append(start_col)
        if end_col is not None:
            usecols.append(end_col)

    if rate_col is not None:
        usecols.append(rate_col)
    else:
        usecols.extend([mod_col, cov_col])

    usecols = list(dict.fromkeys(usecols))

    chunk_bg_files = {}
    total = 0
    kept = 0

    reader = pd.read_csv(infile, sep="\t", usecols=usecols, chunksize=2000000)

    for df in reader:
        total += len(df)

        df = df[df[cov_col] >= args.min_cov].copy()
        if df.empty:
            continue

        if pos_col is not None:
            # 默认pos是1-based，转换成BED 0-based half-open
            df["start0"] = df[pos_col].astype(np.int64) - 1
            df["end0"] = df[pos_col].astype(np.int64)
        else:
            df["start0"] = df[start_col].astype(np.int64)
            if end_col is not None:
                df["end0"] = df[end_col].astype(np.int64)
            else:
                df["end0"] = df["start0"] + 1

        if rate_col is not None:
            df["score"] = df[rate_col].astype(float)
        else:
            df["score"] = df[mod_col].astype(float) / df[cov_col].astype(float)

        df = df[(df["score"] >= 0) & (df["score"] <= 1)]
        if df.empty:
            continue

        df["chunk"] = (df["start0"] // args.chunk_size).astype(np.int64)
        kept += len(df)

        for (chrom, chunk), sub in df.groupby([chrom_col, "chunk"]):
            chunk_start = int(chunk) * args.chunk_size
            chunk_end = chunk_start + args.chunk_size
            bg = outdir / f"{prefix}.{chrom}.{chunk_start}.{chunk_end}.mincov{args.min_cov}.bedGraph"

            sub = sub[[chrom_col, "start0", "end0", "score"]].sort_values([chrom_col, "start0"])
            sub.to_csv(bg, sep="\t", header=False, index=False, mode="a")

            chunk_bg_files[str(bg)] = True

        print(f"[INFO] total={total:,} kept={kept:,}", flush=True)

    print(f"[INFO] generated bedGraph chunks = {len(chunk_bg_files)}", flush=True)

    for bg in sorted(chunk_bg_files):
        bg_path = Path(bg)
        sorted_bg = bg_path.with_suffix(".sorted.bedGraph")
        bw = bg_path.with_suffix(".bw")

        print(f"[INFO] sorting {bg_path.name}", flush=True)
        with open(sorted_bg, "w") as out:
            subprocess.run(
                ["sort", "-k1,1", "-k2,2n", str(bg_path)],
                check=True,
                stdout=out
            )

        print(f"[INFO] converting {sorted_bg.name} -> {bw.name}", flush=True)
        subprocess.run(
            [args.bedGraphToBigWig, str(sorted_bg), args.genome_size, str(bw)],
            check=True
        )

        # 保留sorted bedGraph，删除未排序版本，节约空间
        os.remove(bg_path)

    print("[DONE]", flush=True)

if __name__ == "__main__":
    main()
