#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(".")
IN = BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/P7_marker_module_score/cluster_marker_group_activity_summary.tsv"

OUT = BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/P7_marker_module_score/all_cluster_marker_trend"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(IN, sep="\t")

# 只用包含cluster0-7的版本
d = df[df["version"] == "P7_QCpass_strict"].copy()

if d.empty:
    raise RuntimeError("No P7_QCpass_strict records found in cluster_marker_group_activity_summary.tsv")

# 过滤覆盖太少的module-cluster组合
d = d[
    (d["n_cells"] >= 40) &
    (d["total_cov"] >= 2000)
].copy()

print("[INFO] records after coverage filter:", len(d))
print("[INFO] clusters:", sorted(d["cluster"].unique()))
print("[INFO] groups:", sorted(d["group"].unique()))

# weighted_activity矩阵：group × cluster
mat = d.pivot_table(
    index="group",
    columns="cluster",
    values="weighted_activity",
    aggfunc="mean"
)

# cluster列排序
mat = mat[sorted(mat.columns)]
mat.to_csv(OUT / "P7_QCpass_strict.all_clusters.weighted_activity.matrix.tsv", sep="\t")

# row z-score，只看相对趋势
mat_fill = mat.apply(lambda x: x.fillna(x.mean()), axis=1)
z = mat_fill.sub(mat_fill.mean(axis=1), axis=0)
z = z.div(mat_fill.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
z = z.clip(-2, 2)

z.to_csv(OUT / "P7_QCpass_strict.all_clusters.weighted_activity.row_zscore.tsv", sep="\t")

# 每个module差异幅度
rows = []
for group, row in mat.iterrows():
    s = row.dropna()
    rows.append({
        "group": group,
        "best_cluster": int(s.idxmax()),
        "max_weighted_activity": s.max(),
        "min_weighted_activity": s.min(),
        "delta": s.max() - s.min()
    })

delta = pd.DataFrame(rows).sort_values("delta", ascending=False)
delta.to_csv(OUT / "P7_QCpass_strict.all_clusters.module_delta_summary.tsv", sep="\t", index=False)

# 每个cluster top module
rows = []
for cl in mat.columns:
    s = mat[cl].dropna().sort_values(ascending=False)
    for rank, (group, value) in enumerate(s.head(8).items(), start=1):
        rows.append({
            "cluster": int(cl),
            "rank": rank,
            "group": group,
            "weighted_activity": value,
            "module_delta": float(delta.loc[delta["group"] == group, "delta"].iloc[0])
        })

top = pd.DataFrame(rows)
top.to_csv(OUT / "P7_QCpass_strict.all_clusters.top_modules_by_cluster.tsv", sep="\t", index=False)

# 趋势注释表
anno_rows = []
for cl in mat.columns:
    sub = top[top["cluster"] == cl].sort_values("rank")
    best = sub.iloc[0]
    delta_v = best["module_delta"]

    if delta_v >= 0.05:
        support = "moderate_trend"
    elif delta_v >= 0.02:
        support = "weak_trend"
    elif delta_v >= 0.01:
        support = "very_weak_trend"
    else:
        support = "unclear"

    anno_rows.append({
        "cluster": int(cl),
        "top_module": best["group"],
        "top_weighted_activity": best["weighted_activity"],
        "top_module_delta": delta_v,
        "support": support,
        "top3_modules": ";".join(sub.head(3)["group"].astype(str).tolist())
    })

anno = pd.DataFrame(anno_rows)
anno.to_csv(OUT / "P7_QCpass_strict.all_clusters.trend_annotation.tsv", sep="\t", index=False)

# 热图
plt.figure(figsize=(8, max(4, 0.45 * z.shape[0])))
plt.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
plt.colorbar(label="row z-score of weighted promoter activity")
plt.xticks(range(len(z.columns)), [f"C{x}" for x in z.columns], rotation=45, ha="right")
plt.yticks(range(len(z.index)), z.index)
plt.xlabel("P7 cluster")
plt.ylabel("Pig frontal lobe marker module")
plt.title("P7_QCpass_strict: all cluster marker module trend")
plt.tight_layout()
plt.savefig(OUT / "P7_QCpass_strict.all_clusters.marker_module_heatmap.png", dpi=300)
plt.close()

print("[INFO] wrote:", OUT)
print("\n[INFO] module delta summary:")
print(delta.to_string(index=False))
print("\n[INFO] trend annotation:")
print(anno.to_string(index=False))
