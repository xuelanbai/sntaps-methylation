#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region-counts", required=True)
    ap.add_argument("--cell-list", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunksize", type=int, default=1000000)
    args = ap.parse_args()

    cell_df = pd.read_csv(args.cell_list, sep="\t")
    if "CB" not in cell_df.columns:
        raise ValueError(f"CB column not found in cell-list: {args.cell_list}")

    keep = set(cell_df["CB"].astype(str))
    print(f"[INFO] keep cells: {len(keep)}", flush=True)
    print(f"[INFO] input region-counts: {args.region_counts}", flush=True)
    print(f"[INFO] output: {args.out}", flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    n_in = 0
    n_out = 0
    first = True

    with gzip.open(out, "wt") as fo:
        for chunk in pd.read_csv(args.region_counts, sep="\t", chunksize=args.chunksize):
            if "CB" not in chunk.columns:
                raise ValueError("CB column not found in region-counts")
            chunk["CB"] = chunk["CB"].astype(str)

            n_in += len(chunk)
            chunk = chunk[chunk["CB"].isin(keep)].copy()
            if len(chunk) == 0:
                continue

            n_out += len(chunk)
            chunk.to_csv(fo, sep="\t", index=False, header=first)
            first = False

    print(f"[INFO] input rows: {n_in}", flush=True)
    print(f"[INFO] kept rows: {n_out}", flush=True)

    if n_out == 0:
        raise SystemExit("[ERROR] no rows kept. Check whether CB formats match.")

if __name__ == "__main__":
    main()
