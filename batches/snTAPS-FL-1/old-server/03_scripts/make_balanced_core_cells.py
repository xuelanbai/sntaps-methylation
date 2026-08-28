#!/usr/bin/env python3
import numpy as np
import pandas as pd
from pathlib import Path

qc_path = "08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.global_per_cb.tsv"
outdir = Path("08_qc_summary/final_single_sample/cell_filter_balanced_core")
outdir.mkdir(parents=True, exist_ok=True)

qc = pd.read_csv(qc_path, sep="\t")
print("[INFO] raw CBs:", qc.shape[0])
print("[INFO] columns:", qc.columns.tolist())

for c in ["reads_used", "modified_calls", "covered_calls", "modified_rate", "covered_bins"]:
    qc[c] = pd.to_numeric(qc[c], errors="coerce")

base = qc[
    (qc["reads_used"] >= 100000) &
    (qc["covered_calls"] > 0) &
    (qc["covered_bins"] > 0) &
    (qc["modified_rate"].notna())
].copy()

mlo, mhi = base["modified_rate"].quantile([0.10, 0.90])
base = base[(base["modified_rate"] >= mlo) & (base["modified_rate"] <= mhi)].copy()

# 用robust z-score找“典型细胞”：不是最高覆盖，而是远离低质和极端高覆盖
metrics = pd.DataFrame({
    "log_reads_used": np.log10(base["reads_used"] + 1),
    "log_covered_calls": np.log10(base["covered_calls"] + 1),
    "log_covered_bins": np.log10(base["covered_bins"] + 1),
    "modified_rate": base["modified_rate"]
})

score = np.zeros(base.shape[0], dtype=float)
for c in metrics.columns:
    x = metrics[c].to_numpy()
    med = np.nanmedian(x)
    q25, q75 = np.nanpercentile(x, [25, 75])
    iqr = q75 - q25
    if iqr <= 0:
        iqr = np.nanstd(x)
    if iqr <= 0:
        iqr = 1
    z = (x - med) / iqr
    score += z ** 2

base["balanced_core_score"] = score
base = base.sort_values("balanced_core_score", ascending=True)

print("[INFO] after base filter:", base.shape[0])
print("[INFO] modified_rate q10/q90:", mlo, mhi)

for n in [20000, 10000, 5000, 2000]:
    sub = base.head(n).copy()
    prefix = f"snTAPS-FL-1.balanced_core{len(sub)}"
    sub.to_csv(outdir / f"{prefix}.global_per_cb.tsv", sep="\t", index=False)
    sub["CB"].to_csv(outdir / f"{prefix}.cells.txt", sep="\t", index=False, header=False)
    print("[INFO] wrote", prefix, len(sub))

base.describe().to_csv(outdir / "balanced_core_base_filter.describe.tsv", sep="\t")
