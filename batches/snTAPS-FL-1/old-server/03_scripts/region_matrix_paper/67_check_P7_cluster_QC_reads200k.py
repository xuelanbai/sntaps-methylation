#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path(".")
UMAP = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"
QC = BASE / "08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.global_per_cb.tsv"
OUT = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_reads200k_QC_check"
OUT.mkdir(parents=True, exist_ok=True)

u = pd.read_csv(UMAP, sep="\t")
q = pd.read_csv(QC, sep="\t")

u["CB"] = u["CB"].astype(str)
q["CB"] = q["CB"].astype(str)

qc_cols = ["reads_used", "covered_calls", "covered_bins", "modified_rate"]
keep_qc = ["CB"] + [c for c in qc_cols if c in q.columns]

# 如果UMAP里已有QC列，避免重复；否则从QC表补
for c in qc_cols:
    if c in u.columns:
        u[c] = pd.to_numeric(u[c], errors="coerce")

missing = [c for c in qc_cols if c not in u.columns]
if missing:
    u = u.merge(q[keep_qc], on="CB", how="left", suffixes=("", "_qc"))

for c in qc_cols:
    if c in u.columns:
        u[c] = pd.to_numeric(u[c], errors="coerce")

u["cluster"] = pd.to_numeric(u["cluster"], errors="coerce").astype("Int64")
u["reads_ge200k"] = u["reads_used"] >= 200000

# 每个cluster的QC
summary = (
    u.groupby("cluster")
    .agg(
        n_cells=("CB", "size"),
        reads_ge200k_n=("reads_ge200k", "sum"),
        reads_used_median=("reads_used", "median"),
        reads_used_q10=("reads_used", lambda x: np.nanquantile(x, 0.1)),
        reads_used_q90=("reads_used", lambda x: np.nanquantile(x, 0.9)),
        covered_calls_median=("covered_calls", "median"),
        covered_bins_median=("covered_bins", "median"),
        modified_rate_median=("modified_rate", "median"),
    )
    .reset_index()
)
summary["reads_ge200k_frac"] = summary["reads_ge200k_n"] / summary["n_cells"]
summary.to_csv(OUT / "P7_cluster_QC_summary.tsv", sep="\t", index=False)

# 过滤前后cluster数量
counts_all = u["cluster"].value_counts().sort_index().rename("all_cells")
counts_ge200k = u.loc[u["reads_ge200k"], "cluster"].value_counts().sort_index().rename("reads_ge200k")
counts_lt200k = u.loc[~u["reads_ge200k"], "cluster"].value_counts().sort_index().rename("reads_lt200k")
counts = pd.concat([counts_all, counts_ge200k, counts_lt200k], axis=1).fillna(0).astype(int)
counts.to_csv(OUT / "P7_cluster_counts_before_after_reads200k.tsv", sep="\t")

# 输出过滤后的cell QC，供后续重跑P7使用
q200 = q[q["CB"].isin(set(u.loc[u["reads_ge200k"], "CB"]))].copy()
q200.to_csv(OUT / "snTAPS-FL-1.P7_reads_ge200k.global_per_cb.tsv", sep="\t", index=False)

# 可选：去掉原P7 cluster0/1后的高reads细胞
bad = set(u.loc[u["cluster"].isin([0,1]), "CB"])
q200_no01 = q200[~q200["CB"].isin(bad)].copy()
q200_no01.to_csv(OUT / "snTAPS-FL-1.P7_reads_ge200k_drop_old_cluster0_1.global_per_cb.tsv", sep="\t", index=False)

print("[INFO] wrote:", OUT)
print("[INFO] cells reads>=200k:", len(q200))
print("[INFO] cells reads>=200k and drop old cluster0/1:", len(q200_no01))
