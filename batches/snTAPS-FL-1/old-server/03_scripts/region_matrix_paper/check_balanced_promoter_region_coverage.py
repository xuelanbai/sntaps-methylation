#!/usr/bin/env python3
import gzip
from pathlib import Path
import pandas as pd

files = sorted(Path("08_qc_summary/final_single_sample/region_matrix_balanced_core").glob("*.promoter_TSS_1.5kb.region_counts.tsv.gz"))

for f in files:
    print(f"### {f}")
    usecols = None

    region_counter = {}
    cell_set = set()
    total_rows = 0

    for chunk in pd.read_csv(f, sep="\t", chunksize=1000000):
        total_rows += len(chunk)
        if "CB" not in chunk.columns:
            raise ValueError("CB column not found")
        region_col = "region_id" if "region_id" in chunk.columns else chunk.columns[1]

        chunk["CB"] = chunk["CB"].astype(str)
        chunk[region_col] = chunk[region_col].astype(str)
        cell_set.update(chunk["CB"].unique())

        tmp = chunk.groupby(region_col)["CB"].nunique()
        for r, n in tmp.items():
            region_counter[r] = region_counter.get(r, 0) + int(n)

    s = pd.Series(region_counter)
    print("rows", total_rows)
    print("cells", len(cell_set))
    print("regions_total", len(s))
    for t in [1,3,5,10,20,30,50,100,200,500,1000]:
        print(f"n_cells>={t}", int((s >= t).sum()))
    print()
