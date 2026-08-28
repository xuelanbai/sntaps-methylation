#!/usr/bin/env python3
import argparse
import gzip
import re
from pathlib import Path
import pandas as pd
import numpy as np

def open_text(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def parse_bin_id(x):
    s = str(x)
    # 常见格式：1:100000-200000
    m = re.match(r"^([^:\s]+):(\d+)-(\d+)$", s)
    if m:
        chrom = m.group(1)
        start = int(m.group(2))
        end = int(m.group(3))
        s500 = (start // 500000) * 500000
        e500 = s500 + 500000
        return f"{chrom}:{s500}-{e500}"

    # 如果是chrom_start_end
    m = re.match(r"^([^_\s]+)_(\d+)_(\d+)$", s)
    if m:
        chrom = m.group(1)
        start = int(m.group(2))
        s500 = (start // 500000) * 500000
        e500 = s500 + 500000
        return f"{chrom}:{s500}-{e500}"

    raise ValueError(f"Cannot parse bin id: {s}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin-counts", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--feature-bin-list", default=None)
    ap.add_argument("--chunksize", type=int, default=1000000)
    args = ap.parse_args()

    allowed = None
    if args.feature_bin_list:
        allowed = set()
        with open(args.feature_bin_list) as f:
            for line in f:
                x = line.strip()
                if x:
                    allowed.add(x)
        print(f"[INFO] feature-bin-list: {args.feature_bin_list}, n={len(allowed)}", flush=True)

    print("[INFO] inspect header", flush=True)
    with open_text(args.bin_counts) as f:
        header = f.readline().rstrip("\n").split("\t")
    print("[INFO] columns:", header, flush=True)

    cb_col = "CB" if "CB" in header else header[0]
    bin_col = None
    for c in ["bin_id", "region_id", "bin", "feature_id"]:
        if c in header:
            bin_col = c
            break
    if bin_col is None:
        # 常见情况第二列是bin
        bin_col = header[1]

    mod_col = "modified_calls"
    cov_col = "covered_calls"
    if mod_col not in header or cov_col not in header:
        raise ValueError("modified_calls/covered_calls not found in bin_counts")

    tmp = Path(str(args.out) + ".tmp.tsv")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    print("[INFO] first pass: chunk aggregate to tmp", flush=True)
    first = True
    usecols = [cb_col, bin_col, mod_col, cov_col]
    for chunk in pd.read_csv(args.bin_counts, sep="\t", usecols=usecols, chunksize=args.chunksize):
        chunk[cb_col] = chunk[cb_col].astype(str)
        chunk[bin_col] = chunk[bin_col].astype(str)

        if allowed is not None:
            chunk = chunk[chunk[bin_col].isin(allowed)].copy()
            if chunk.empty:
                continue

        chunk[mod_col] = pd.to_numeric(chunk[mod_col], errors="coerce").fillna(0)
        chunk[cov_col] = pd.to_numeric(chunk[cov_col], errors="coerce").fillna(0)
        chunk = chunk[chunk[cov_col] > 0].copy()
        if chunk.empty:
            continue

        chunk["bin500_id"] = chunk[bin_col].map(parse_bin_id)
        g = chunk.groupby([cb_col, "bin500_id"], as_index=False)[[mod_col, cov_col]].sum()
        g.to_csv(tmp, sep="\t", index=False, header=first, mode="w" if first else "a")
        first = False

    if first:
        raise ValueError("No rows written to tmp; check feature list/input")

    print("[INFO] second pass: final aggregate", flush=True)
    acc = {}
    for chunk in pd.read_csv(tmp, sep="\t", chunksize=args.chunksize):
        chunk[cb_col] = chunk[cb_col].astype(str)
        chunk["bin500_id"] = chunk["bin500_id"].astype(str)
        g = chunk.groupby([cb_col, "bin500_id"], as_index=False)[[mod_col, cov_col]].sum()
        for cb, bid, m, c in g[[cb_col, "bin500_id", mod_col, cov_col]].itertuples(index=False):
            key = (cb, bid)
            if key not in acc:
                acc[key] = [0.0, 0.0]
            acc[key][0] += float(m)
            acc[key][1] += float(c)

    print("[INFO] write output", out, flush=True)
    with gzip.open(out, "wt") as fo:
        fo.write("CB\tbin_id\tmodified_calls\tcovered_calls\tmodified_rate\n")
        for (cb, bid), (m, c) in acc.items():
            if c <= 0:
                continue
            fo.write(f"{cb}\t{bid}\t{int(round(m))}\t{int(round(c))}\t{m/c:.8f}\n")

    tmp.unlink(missing_ok=True)
    print("[INFO] done", flush=True)

if __name__ == "__main__":
    main()
