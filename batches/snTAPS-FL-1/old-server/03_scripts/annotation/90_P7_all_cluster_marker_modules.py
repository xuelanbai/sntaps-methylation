#!/usr/bin/env python3
from pathlib import Path
import gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(".")

UMAP = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"

REGION_MATRIX = BASE / "08_qc_summary/final_single_sample/region_matrix_covfilter/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.cov10.region_counts.tsv.gz"

OUTDIR = BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/P7_all_cluster_marker_module"
OUTDIR.mkdir(parents=True, exist_ok=True)

# 优先找之前已经生成过的猪额叶marker-promoter对应文件
candidate_marker_files = [
    BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/pig_frontal_lobe_marker_promoters.tsv",
    BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/P7_marker_module_score/pig_frontal_lobe_marker_promoters.tsv",
    BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/P7_marker_module_score/marker_promoters.tsv",
    BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/marker_promoters.tsv",
]

marker_file = None
for f in candidate_marker_files:
    if f.exists():
        marker_file = f
        break

if marker_file is None:
    hits = list((BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe").rglob("*marker*promoter*.tsv"))
    if len(hits) > 0:
        marker_file = hits[0]

if marker_file is None:
    raise FileNotFoundError("Cannot find marker-promoter mapping tsv under marker_pig_frontal_lobe")

print("[INFO] marker file:", marker_file, flush=True)

umap = pd.read_csv(UMAP, sep="\t")
print("[INFO] UMAP columns:", umap.columns.tolist(), flush=True)

# 自动识别CB和cluster列
cb_col = "CB" if "CB" in umap.columns else ("cell" if "cell" in umap.columns else umap.columns[0])
cluster_col = None
for c in ["cluster", "leiden", "kmeans", "P7_cluster"]:
    if c in umap.columns:
        cluster_col = c
        break
if cluster_col is None:
    # 常见情况：最后一列是cluster
    cluster_col = umap.columns[-1]

umap = umap[[cb_col, cluster_col]].rename(columns={cb_col: "CB", cluster_col: "cluster"})
umap["cluster"] = umap["cluster"].astype(str)

print("[INFO] cluster counts:")
print(umap["cluster"].value_counts().sort_index().to_string(), flush=True)

markers = pd.read_csv(marker_file, sep="\t")
print("[INFO] marker columns:", markers.columns.tolist(), flush=True)

# 兼容不同列名
rename = {}
for c in markers.columns:
    lc = c.lower()
    if lc in ["group", "celltype", "cell_type", "type", "module"]:
        rename[c] = "group"
    if lc in ["gene", "gene_name", "symbol"]:
        rename[c] = "gene"
    if lc in ["region_id", "promoter", "promoter_id"]:
        rename[c] = "region_id"

markers = markers.rename(columns=rename)

if "group" not in markers.columns:
    raise ValueError(f"marker file lacks group column: {markers.columns.tolist()}")
if "region_id" not in markers.columns:
    raise ValueError(f"marker file lacks region_id column: {markers.columns.tolist()}")

markers = markers[["group", "region_id"]].dropna().drop_duplicates()
marker_regions = set(markers["region_id"].astype(str))
print("[INFO] marker groups:", markers["group"].nunique(), "marker regions:", len(marker_regions), flush=True)

# 准备三套细胞集合：
# 1 all：C0-C7全部
# 2 no01_reference：只看C2-C7
# 3 no016_reference：去掉C0/C1/C6
sets = {
    "P7_all_clusters": set(umap["CB"]),
    "P7_no01_reference": set(umap.loc[~umap["cluster"].isin(["0", "1"]), "CB"]),
    "P7_no016_reference": set(umap.loc[~umap["cluster"].isin(["0", "1", "6"]), "CB"]),
}

# 如果已有QCpass cell list，也自动加入
qcpass_candidates = [
    BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_QCpass_old_umap/P7_QCpass_strict.cells.txt",
    BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_remove_abnormal_old_umap/P7_info_strict.cells.txt",
]
for f in qcpass_candidates:
    if f.exists():
        cells = set(pd.read_csv(f, header=None)[0].astype(str))
        sets[f.stem] = cells
        print("[INFO] add cell set from", f, len(cells), flush=True)

def read_marker_rows(matrix_path, selected_cells):
    use_rows = []
    with gzip.open(matrix_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {x:i for i,x in enumerate(header)}
        required = ["CB", "region_id", "modified_calls", "covered_calls"]
        for r in required:
            if r not in idx:
                raise ValueError(f"{r} not in matrix header: {header}")
        for line in fh:
            arr = line.rstrip("\n").split("\t")
            cb = arr[idx["CB"]]
            if cb not in selected_cells:
                continue
            rid = arr[idx["region_id"]]
            if rid not in marker_regions:
                continue
            try:
                mod = float(arr[idx["modified_calls"]])
                cov = float(arr[idx["covered_calls"]])
            except Exception:
                continue
            if cov <= 0:
                continue
            use_rows.append((cb, rid, mod, cov))
    return pd.DataFrame(use_rows, columns=["CB", "region_id", "modified_calls", "covered_calls"])

def zscore_rows(mat):
    x = mat.to_numpy(dtype=float)
    mu = np.nanmean(x, axis=1, keepdims=True)
    sd = np.nanstd(x, axis=1, keepdims=True)
    sd[sd == 0] = np.nan
    z = (x - mu) / sd
    z = np.nan_to_num(z, nan=0.0)
    z = np.clip(z, -2, 2)
    return pd.DataFrame(z, index=mat.index, columns=mat.columns)

for name, cells in sets.items():
    print(f"[INFO] process {name}: selected cells={len(cells)}", flush=True)
    sub_umap = umap[umap["CB"].isin(cells)].copy()
    rows = read_marker_rows(REGION_MATRIX, cells)
    print(f"[INFO] marker matrix rows={len(rows)}", flush=True)

    if rows.empty:
        print(f"[WARN] no marker rows for {name}", flush=True)
        continue

    rows = rows.merge(sub_umap, on="CB", how="inner")
    rows = rows.merge(markers, on="region_id", how="inner")

    # 加权promoter低甲基化活性：1 - sum(modified_calls)/sum(covered_calls)
    g = rows.groupby(["group", "cluster"], as_index=False).agg(
        modified_calls=("modified_calls", "sum"),
        covered_calls=("covered_calls", "sum"),
        n_obs=("CB", "count"),
        n_cells=("CB", "nunique"),
        n_regions=("region_id", "nunique"),
    )
    g["weighted_modified_rate"] = g["modified_calls"] / g["covered_calls"]
    g["weighted_activity"] = 1 - g["weighted_modified_rate"]

    g.to_csv(OUTDIR / f"{name}.module_activity_long.tsv", sep="\t", index=False)

    mat = g.pivot(index="group", columns="cluster", values="weighted_activity")
    mat = mat.reindex(columns=sorted(mat.columns, key=lambda x: int(x) if str(x).isdigit() else str(x)))
    mat.to_csv(OUTDIR / f"{name}.mean_activity.matrix.tsv", sep="\t")

    zmat = zscore_rows(mat)
    zmat.to_csv(OUTDIR / f"{name}.row_zscore.matrix.tsv", sep="\t")

    # 每个cluster最高的module
    top = g.sort_values(["cluster", "weighted_activity"], ascending=[True, False]).groupby("cluster").head(5)
    top.to_csv(OUTDIR / f"{name}.top_modules_by_cluster.tsv", sep="\t", index=False)

    # 每个module的cluster差异
    delta_rows = []
    for grp, df in g.groupby("group"):
        df = df.sort_values("weighted_activity", ascending=False)
        delta = df["weighted_activity"].max() - df["weighted_activity"].min()
        delta_rows.append({
            "group": grp,
            "best_cluster": df.iloc[0]["cluster"],
            "best_activity": df.iloc[0]["weighted_activity"],
            "worst_cluster": df.iloc[-1]["cluster"],
            "worst_activity": df.iloc[-1]["weighted_activity"],
            "delta": delta,
        })
    delta = pd.DataFrame(delta_rows).sort_values("delta", ascending=False)
    delta.to_csv(OUTDIR / f"{name}.module_delta_summary.tsv", sep="\t", index=False)

    # 画热图
    fig_h = max(6, 0.25 * len(zmat))
    plt.figure(figsize=(8, fig_h))
    im = plt.imshow(zmat.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    plt.colorbar(im, label="row z-score of promoter activity")
    plt.xticks(range(len(zmat.columns)), [f"C{x}" for x in zmat.columns], rotation=45)
    plt.yticks(range(len(zmat.index)), zmat.index)
    plt.title(f"{name}: pig frontal marker promoter activity")
    plt.xlabel("P7 cluster")
    plt.ylabel("marker module")
    plt.tight_layout()
    plt.savefig(OUTDIR / f"{name}.marker_module_heatmap.png", dpi=200)
    plt.close()

    print("[INFO] wrote", name, flush=True)

print("[INFO] done:", OUTDIR)
