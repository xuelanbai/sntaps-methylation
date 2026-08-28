#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(".")
INDIR = BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/P7_marker_module_score"
IN = INDIR / "cluster_marker_group_activity_summary.tsv"
OUT = INDIR / "weighted_replot"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(IN, sep="\t")

# 只看覆盖足够的module-cluster组合，太少的组合不用于判断
MIN_N_CELLS = 40
MIN_TOTAL_COV = 2000

df["pass_cov"] = (df["n_cells"] >= MIN_N_CELLS) & (df["total_cov"] >= MIN_TOTAL_COV)

for version, d0 in df.groupby("version"):
    d = d0[d0["pass_cov"]].copy()

    for value_col in ["weighted_activity", "mean_activity"]:
        mat = d.pivot_table(index="group", columns="cluster", values=value_col, aggfunc="mean")
        nmat = d.pivot_table(index="group", columns="cluster", values="n_cells", aggfunc="sum")

        mat.to_csv(OUT / f"{version}.{value_col}.matrix.tsv", sep="\t")
        nmat.to_csv(OUT / f"{version}.{value_col}.n_cells.tsv", sep="\t")

        # row z-score，仅用于看相对趋势
        mat_fill = mat.apply(lambda x: x.fillna(x.mean()), axis=1)
        z = mat_fill.sub(mat_fill.mean(axis=1), axis=0)
        z = z.div(mat_fill.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
        z.to_csv(OUT / f"{version}.{value_col}.row_zscore.tsv", sep="\t")

        plt.figure(figsize=(7, max(4, 0.45 * z.shape[0])))
        im = plt.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
        plt.xticks(range(z.shape[1]), [f"C{x}" for x in z.columns], rotation=45, ha="right")
        plt.yticks(range(z.shape[0]), z.index)
        plt.xlabel("P7 cluster")
        plt.ylabel("Marker module")
        plt.title(f"{version}: {value_col} row z-score")
        plt.colorbar(im, label="row z-score")
        plt.tight_layout()
        plt.savefig(OUT / f"{version}.{value_col}.row_zscore_heatmap.png", dpi=240)
        plt.close()

    # 输出每个cluster最高的几个module，只作为探索
    w = d.pivot_table(index="group", columns="cluster", values="weighted_activity", aggfunc="mean")
    rows = []
    for cl in w.columns:
        s = w[cl].dropna().sort_values(ascending=False)
        for rank, (group, score) in enumerate(s.head(5).items(), start=1):
            rows.append({
                "cluster": cl,
                "rank": rank,
                "group": group,
                "weighted_activity": score
            })
    pd.DataFrame(rows).to_csv(OUT / f"{version}.top_modules_by_cluster.tsv", sep="\t", index=False)

    # 输出每个module在哪个cluster最高，以及最大-最小差值
    rows = []
    for group, s in w.iterrows():
        s = s.dropna()
        if len(s) == 0:
            continue
        rows.append({
            "group": group,
            "best_cluster": s.idxmax(),
            "max_weighted_activity": s.max(),
            "min_weighted_activity": s.min(),
            "delta": s.max() - s.min()
        })
    pd.DataFrame(rows).sort_values("delta", ascending=False).to_csv(
        OUT / f"{version}.module_delta_summary.tsv", sep="\t", index=False
    )

print("[INFO] done")
