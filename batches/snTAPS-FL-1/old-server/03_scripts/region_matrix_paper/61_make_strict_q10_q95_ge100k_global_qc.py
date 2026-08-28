#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

base = Path("08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.global_per_cb.tsv")
outdir = Path("08_qc_summary/final_single_sample/cell_filter_strict_q10_q95_ge100k")
outdir.mkdir(parents=True, exist_ok=True)

out = outdir / "snTAPS-FL-1.strict_q10_q95_ge100k.global_per_cb.tsv"
summary_out = outdir / "snTAPS-FL-1.strict_q10_q95_ge100k.filter_summary.tsv"
cb_out = outdir / "snTAPS-FL-1.strict_q10_q95_ge100k.CBlist.txt"

if not base.exists():
    raise FileNotFoundError(base)

qc = pd.read_csv(base, sep="\t")
print("[INFO] input QC:", base)
print("[INFO] columns:", list(qc.columns))
print("[INFO] input cells:", qc.shape[0])

# 自动识别CB列
cb_candidates = ["CB", "cell", "barcode", "raw_CB"]
cb_col = next((c for c in cb_candidates if c in qc.columns), None)
if cb_col is None:
    raise ValueError(f"Cannot find CB column. columns={list(qc.columns)}")

# 自动识别reads列
read_candidates = [
    "reads_used", "unique_reads", "mapped_unique_reads",
    "n_unique_reads", "total_reads", "reads", "read_count",
    "n_reads", "num_reads", "final_reads", "dedup_reads"
]
read_col = next((c for c in read_candidates if c in qc.columns), None)

if read_col is None:
    raise ValueError(
        "Cannot find reads column in global QC. "
        "Please check columns above. Need one of: " + ",".join(read_candidates)
    )

# 统一成reads_used，方便后续脚本使用
if read_col != "reads_used":
    qc["reads_used"] = pd.to_numeric(qc[read_col], errors="coerce")
else:
    qc["reads_used"] = pd.to_numeric(qc["reads_used"], errors="coerce")

qc = qc.dropna(subset=["reads_used"]).copy()

q10 = qc["reads_used"].quantile(0.10)
q95 = qc["reads_used"].quantile(0.95)
lower = max(q10, 100000)

keep = (qc["reads_used"] >= lower) & (qc["reads_used"] <= q95)
filt = qc.loc[keep].copy()

# 保证CB列名叫CB
if cb_col != "CB":
    filt = filt.rename(columns={cb_col: "CB"})

filt.to_csv(out, sep="\t", index=False)
filt["CB"].drop_duplicates().to_csv(cb_out, sep="\t", index=False, header=False)

summary = pd.DataFrame([{
    "input_cells": qc.shape[0],
    "read_column_used": read_col,
    "q10_reads": q10,
    "q95_reads": q95,
    "lower_cutoff_used": lower,
    "upper_cutoff_used": q95,
    "output_cells": filt.shape[0],
}])
summary.to_csv(summary_out, sep="\t", index=False)

print("[INFO] read column used:", read_col)
print("[INFO] q10:", q10)
print("[INFO] q95:", q95)
print("[INFO] lower cutoff used:", lower)
print("[INFO] output cells:", filt.shape[0])
print("[INFO] output:", out)
print("[INFO] CBlist:", cb_out)
print("[INFO] summary:", summary_out)
