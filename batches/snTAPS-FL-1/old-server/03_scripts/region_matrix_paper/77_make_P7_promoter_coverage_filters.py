#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path(".")
P7 = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"

PROM = BASE / "08_qc_summary/final_single_sample/region_matrix/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.region_counts.tsv.gz"

OUT = BASE / "08_qc_summary/final_single_sample/cell_filter_P7_promoter_coverage"
OUT.mkdir(parents=True, exist_ok=True)

p7 = pd.read_csv(P7, sep="\t")
p7["CB"] = p7["CB"].astype(str)
keep = set(p7["CB"])

print("[INFO] P7 cells:", len(p7), flush=True)
print("[INFO] scan promoter matrix:", PROM, flush=True)

agg = {}

for chunk in pd.read_csv(PROM, sep="\t", chunksize=1000000):
    chunk["CB"] = chunk["CB"].astype(str)
    chunk = chunk[chunk["CB"].isin(keep)].copy()
    if chunk.empty:
        continue

    chunk["covered_calls"] = pd.to_numeric(chunk["covered_calls"], errors="coerce").fillna(0)
    chunk["modified_calls"] = pd.to_numeric(chunk["modified_calls"], errors="coerce").fillna(0)

    g = chunk.groupby("CB").agg(
        promoter_n_regions=("region_id", "nunique"),
        promoter_covered_calls=("covered_calls", "sum"),
        promoter_modified_calls=("modified_calls", "sum"),
    )

    for cb, row in g.iterrows():
        if cb not in agg:
            agg[cb] = [0, 0.0, 0.0]
        agg[cb][0] += row["promoter_n_regions"]
        agg[cb][1] += row["promoter_covered_calls"]
        agg[cb][2] += row["promoter_modified_calls"]

prom = pd.DataFrame.from_dict(
    agg,
    orient="index",
    columns=["promoter_n_regions", "promoter_covered_calls", "promoter_modified_calls"]
).reset_index().rename(columns={"index": "CB"})

prom["CB"] = prom["CB"].astype(str)
prom["promoter_modified_rate"] = prom["promoter_modified_calls"] / prom["promoter_covered_calls"].replace(0, np.nan)

df = p7.merge(prom, on="CB", how="left")
for c in ["promoter_n_regions", "promoter_covered_calls", "promoter_modified_calls"]:
    df[c] = df[c].fillna(0)
df["promoter_modified_rate"] = df["promoter_modified_rate"].fillna(0)

qc_cols = ["reads_used", "covered_calls", "covered_bins", "modified_rate"]
prom_cols = ["promoter_n_regions", "promoter_covered_calls"]

for c in qc_cols + prom_cols + ["promoter_modified_rate"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df.to_csv(OUT / "snTAPS-FL-1.P7_with_promoter_coverage_qc.tsv", sep="\t", index=False)

def write_filter(name, mask):
    sub = df[mask].copy()
    sub.to_csv(OUT / f"snTAPS-FL-1.{name}.global_per_cb.tsv", sep="\t", index=False)

    row = {"filter": name, "n_cells": len(sub)}
    if "cluster" in sub.columns:
        row["cluster0"] = int((sub["cluster"] == 0).sum())
        row["cluster1"] = int((sub["cluster"] == 1).sum())
        row["cluster0_1"] = row["cluster0"] + row["cluster1"]
        row["cluster0_1_frac"] = row["cluster0_1"] / len(sub) if len(sub) else np.nan

    for c in qc_cols + prom_cols + ["promoter_modified_rate"]:
        row[f"{c}_median"] = sub[c].median()
        row[f"{c}_min"] = sub[c].min()
        row[f"{c}_max"] = sub[c].max()
    return row

summary = []

# A: 只去掉高覆盖异常，保守
m = pd.Series(True, index=df.index)
for c in ["reads_used", "covered_calls", "covered_bins", "promoter_n_regions", "promoter_covered_calls"]:
    m &= df[c] <= df[c].quantile(0.95)
m &= df["modified_rate"].between(df["modified_rate"].quantile(0.05), df["modified_rate"].quantile(0.95))
m &= df["promoter_n_regions"] > 0
summary.append(write_filter("P7_noHigh95", m))

# B: promoter覆盖和global QC都保留中间80%
m = pd.Series(True, index=df.index)
for c in ["reads_used", "covered_calls", "covered_bins", "modified_rate", "promoter_n_regions", "promoter_covered_calls"]:
    m &= df[c].between(df[c].quantile(0.10), df[c].quantile(0.90))
m &= df["promoter_n_regions"] > 0
summary.append(write_filter("P7_covCore80", m))

# C: promoter覆盖保留中间50%，最严格
m = pd.Series(True, index=df.index)
for c in ["reads_used", "covered_calls", "covered_bins", "modified_rate", "promoter_n_regions", "promoter_covered_calls"]:
    m &= df[c].between(df[c].quantile(0.25), df[c].quantile(0.75))
m &= df["promoter_n_regions"] > 0
summary.append(write_filter("P7_covCore50", m))

pd.DataFrame(summary).to_csv(OUT / "P7_promoter_coverage_filter_summary.tsv", sep="\t", index=False)

print("[INFO] summary")
print(pd.DataFrame(summary), flush=True)
