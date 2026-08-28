#!/usr/bin/env python3
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A_UMAP = Path(os.environ.get(
    "A_UMAP",
    "08_qc_summary/final_single_sample/region_umap_paper/P7_no01_promoter_residual_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"
))

REGION_MAP = Path("08_qc_summary/final_single_sample/marker_pig_frontal_lobe/P7_QCpass_marker/marker_promoter_region_map.tsv")
MATRIX = Path("08_qc_summary/final_single_sample/region_matrix/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.region_counts.tsv.gz")

OUTDIR = Path("08_qc_summary/final_single_sample/marker_pig_frontal_lobe/A_no01_rerun_marker_existingMap")
OUTDIR.mkdir(parents=True, exist_ok=True)

print("[INFO] A_UMAP:", A_UMAP, flush=True)
print("[INFO] REGION_MAP:", REGION_MAP, flush=True)
print("[INFO] MATRIX:", MATRIX, flush=True)

umap = pd.read_csv(A_UMAP, sep="\t")
umap["CB"] = umap["CB"].astype(str)
umap["cluster"] = pd.to_numeric(umap["cluster"], errors="coerce").astype(int)

cell2cluster = dict(zip(umap["CB"], umap["cluster"]))
cell_set = set(cell2cluster)

print("[INFO] selected cells:", len(cell_set), flush=True)
print("[INFO] cluster counts:", flush=True)
print(umap["cluster"].value_counts().sort_index().to_string(), flush=True)

region_map = pd.read_csv(REGION_MAP, sep="\t")
region_map["region_id"] = region_map["region_id"].astype(str)

# 兼容不同列名
if "label" not in region_map.columns:
    if {"group", "gene"}.issubset(region_map.columns):
        region_map["label"] = region_map["group"].astype(str) + ":" + region_map["gene"].astype(str)
    else:
        raise ValueError(f"REGION_MAP columns={region_map.columns.tolist()}, need group/gene or label")

if "group" not in region_map.columns:
    region_map["group"] = region_map["label"].astype(str).str.split(":", n=1).str[0]
if "gene" not in region_map.columns:
    region_map["gene"] = region_map["label"].astype(str).str.split(":", n=1).str[1]

region_map = region_map[["region_id", "group", "gene", "label"]].drop_duplicates()
region_map.to_csv(OUTDIR / "A_no01_marker_region_map_existingMap.tsv", sep="\t", index=False)

region2labels = region_map.groupby("region_id")["label"].apply(list).to_dict()
region_set = set(region2labels)

print("[INFO] marker regions:", len(region_set), flush=True)
print("[INFO] marker genes:", region_map["gene"].nunique(), flush=True)

agg = {}
usecols = ["CB", "region_id", "modified_calls", "covered_calls"]

total_cb_hit = 0
total_region_hit = 0
total_both = 0

for i, chunk in enumerate(pd.read_csv(MATRIX, sep="\t", usecols=usecols, chunksize=1000000), 1):
    chunk["CB"] = chunk["CB"].astype(str)
    chunk["region_id"] = chunk["region_id"].astype(str)

    c = chunk["CB"].isin(cell_set)
    r = chunk["region_id"].isin(region_set)

    total_cb_hit += int(c.sum())
    total_region_hit += int(r.sum())

    sub = chunk[c & r].copy()
    total_both += len(sub)

    if sub.empty:
        if i % 10 == 0:
            print(f"[INFO] chunk={i}, total_both={total_both}", flush=True)
        continue

    sub["cluster"] = sub["CB"].map(cell2cluster)
    sub["modified_calls"] = pd.to_numeric(sub["modified_calls"], errors="coerce").fillna(0)
    sub["covered_calls"] = pd.to_numeric(sub["covered_calls"], errors="coerce").fillna(0)

    rows = []
    for _, x in sub.iterrows():
        for lab in region2labels.get(x["region_id"], []):
            rows.append((int(x["cluster"]), lab, float(x["modified_calls"]), float(x["covered_calls"])))

    tmp = pd.DataFrame(rows, columns=["cluster", "label", "modified_calls", "covered_calls"])
    g = tmp.groupby(["cluster", "label"], as_index=False)[["modified_calls", "covered_calls"]].sum()

    for _, rr in g.iterrows():
        key = (int(rr["cluster"]), rr["label"])
        agg.setdefault(key, [0.0, 0.0])
        agg[key][0] += float(rr["modified_calls"])
        agg[key][1] += float(rr["covered_calls"])

    print(f"[INFO] chunk={i}, kept={len(sub)}, total_both={total_both}", flush=True)

print("[INFO] rows_with_A_cells:", total_cb_hit, flush=True)
print("[INFO] rows_with_marker_regions:", total_region_hit, flush=True)
print("[INFO] rows_with_A_cells_AND_marker_regions:", total_both, flush=True)

rows = []
for (cluster, label), (mod, cov) in agg.items():
    if cov <= 0:
        continue
    group, gene = label.split(":", 1)
    rate = mod / cov
    rows.append({
        "cluster": cluster,
        "label": label,
        "group": group,
        "gene": gene,
        "modified_calls": mod,
        "covered_calls": cov,
        "modified_rate": rate,
        "hypomethylation_activity": -rate
    })

long = pd.DataFrame(rows)
long.to_csv(OUTDIR / "A_no01_marker_promoter_activity_long.tsv", sep="\t", index=False)

if long.empty:
    raise RuntimeError("Still no marker records found. Check existing REGION_MAP compatibility.")

# group-level aggregation
group_long = long.groupby(["cluster", "group"], as_index=False)[["modified_calls", "covered_calls"]].sum()
group_long["modified_rate"] = group_long["modified_calls"] / group_long["covered_calls"]
group_long["hypomethylation_activity"] = -group_long["modified_rate"]
group_long.to_csv(OUTDIR / "A_no01_marker_group_activity_long.tsv", sep="\t", index=False)

cluster_order = sorted(group_long["cluster"].unique())
group_order = [
    "Excitatory_neuron", "Upper_layer_excitatory", "Deep_layer_excitatory",
    "Inhibitory_neuron", "Astrocyte", "Oligodendrocyte", "OPC",
    "Microglia", "Endothelial", "Pericyte_SMC"
]

def row_zscore(mat):
    arr = mat.to_numpy(dtype=float)
    mu = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, keepdims=True)
    sd[sd == 0] = np.nan
    z = (arr - mu) / sd
    z = np.nan_to_num(z, nan=0)
    z = np.clip(z, -2, 2)
    return pd.DataFrame(z, index=mat.index, columns=[f"C{c}" for c in mat.columns])

gmat = group_long.pivot_table(index="group", columns="cluster", values="hypomethylation_activity", aggfunc="mean")
gmat = gmat.reindex(index=[g for g in group_order if g in gmat.index], columns=cluster_order)

gz = row_zscore(gmat)
gz.to_csv(OUTDIR / "A_no01_marker_group_activity_zscore.tsv", sep="\t")

# gene-level heatmap
mat = long.pivot_table(index="label", columns="cluster", values="hypomethylation_activity", aggfunc="mean")
mat = mat.reindex(columns=cluster_order)
z = row_zscore(mat)
z.to_csv(OUTDIR / "A_no01_marker_promoter_activity_zscore.tsv", sep="\t")

# delta summary
delta_rows = []
for group, row in gmat.iterrows():
    s = row.dropna()
    delta_rows.append({
        "group": group,
        "best_cluster": int(s.idxmax()),
        "max_activity": s.max(),
        "min_activity": s.min(),
        "delta": s.max() - s.min()
    })

delta = pd.DataFrame(delta_rows).sort_values("delta", ascending=False)
delta.to_csv(OUTDIR / "A_no01_marker_group_delta_summary.tsv", sep="\t", index=False)

def plot_heatmap(df, outprefix, title, figsize):
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(df.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(df.shape[1]))
    ax.set_xticklabels(df.columns)
    ax.set_yticks(range(df.shape[0]))
    ax.set_yticklabels(df.index, fontsize=7 if df.shape[0] > 20 else 10)
    ax.set_xlabel("A rerun cluster")
    ax.set_title(title)
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Relative promoter hypomethylation activity\nrow z-score")
    plt.tight_layout()
    fig.savefig(str(outprefix) + ".png", dpi=300)
    fig.savefig(str(outprefix) + ".pdf")
    plt.close(fig)

plot_heatmap(
    gz,
    OUTDIR / "A_no01_marker_group_activity_heatmap",
    "A rerun UMAP: marker group activity",
    (7.2, 4.8)
)

plot_heatmap(
    z,
    OUTDIR / "A_no01_marker_promoter_activity_heatmap",
    "A rerun UMAP: promoter hypomethylation marker activity",
    (8, max(8, 0.23 * z.shape[0]))
)

top_rows = []
for c in gz.columns:
    s = gz[c].sort_values(ascending=False)
    for rank, (grp, score) in enumerate(s.head(5).items(), 1):
        top_rows.append({
            "cluster": c,
            "rank": rank,
            "marker_group": grp,
            "zscore": score
        })

top = pd.DataFrame(top_rows)
top.to_csv(OUTDIR / "A_no01_top_marker_groups_by_cluster.tsv", sep="\t", index=False)

print("[INFO] done:", OUTDIR, flush=True)
print("[INFO] top marker groups:", flush=True)
print(top.query("rank <= 3").to_string(index=False), flush=True)
print("[INFO] delta summary:", flush=True)
print(delta.to_string(index=False), flush=True)
