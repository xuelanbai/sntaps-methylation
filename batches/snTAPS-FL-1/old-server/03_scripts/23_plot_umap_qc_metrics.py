#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umap", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.umap, sep="\t")

    xcol = "UMAP1" if "UMAP1" in df.columns else "PC1"
    ycol = "UMAP2" if "UMAP2" in df.columns else "PC2"

    metrics = ["reads_used", "covered_calls", "covered_bins", "modified_rate"]

    for m in metrics:
        if m not in df.columns:
            continue

        plt.figure(figsize=(6,5))
        sc = plt.scatter(df[xcol], df[ycol], c=df[m], s=6, alpha=0.85)
        plt.xlabel(xcol)
        plt.ylabel(ycol)
        plt.title(f"{args.prefix} colored by {m}")
        cbar = plt.colorbar(sc)
        cbar.set_label(m)
        plt.tight_layout()
        plt.savefig(outdir / f"{args.prefix}.{m}.png", dpi=300)
        plt.close()

    print("[INFO] done")

if __name__ == "__main__":
    main()
