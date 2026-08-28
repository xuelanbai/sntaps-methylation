#!/usr/bin/env python3
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(".")
MARKER_TSV = BASE / "03_scripts/annotation/pig_frontal_lobe_marker_symbols.tsv"
BED = BASE / "02_reference/pig/gene_regions/pig_promoter_TSS_1.5kb.bed"
REGION_COUNTS = BASE / "08_qc_summary/final_single_sample/region_matrix_covfilter/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.cov10.region_counts.tsv.gz"

P7_UMAP = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"

CELL_FILES = {
    "P7_QCpass_strict": BASE / "08_qc_summary/final_single_sample/cell_filter_P7_QCpass/snTAPS-FL-1.P7_QCpass_strict.global_per_cb.tsv",
    "P7_QCpass_no01_reference": BASE / "08_qc_summary/final_single_sample/cell_filter_P7_QCpass/snTAPS-FL-1.P7_QCpass_no01_reference.global_per_cb.tsv",
    "P7_QCpass_no016_reference": BASE / "08_qc_summary/final_single_sample/cell_filter_P7_QCpass/snTAPS-FL-1.P7_QCpass_no016_reference.global_per_cb.tsv",
}

OUT = BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/P7_marker_module_score"
OUT.mkdir(parents=True, exist_ok=True)

MIN_GROUP_COV_PER_CELL = 30

print("[INFO] read markers", flush=True)
markers = pd.read_csv(MARKER_TSV, sep="\t")
markers["gene"] = markers["gene"].astype(str).str.strip()
markers["group"] = markers["group"].astype(str).str.strip()
markers["gene_upper"] = markers["gene"].str.upper()
marker_set = set(markers["gene_upper"])

print("[INFO] parse promoter BED", flush=True)
rows = []
with open(BED) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = re.split(r"\s+", line.strip())
        ens = [x for x in parts if re.match(r"^ENSSSCG\d+", x)]
        if not ens:
            continue
        region_id = ens[0]
        for x in parts:
            xu = x.upper()
            if xu in marker_set:
                rows.append((region_id, xu))

map_df = pd.DataFrame(rows, columns=["region_id", "gene_upper"]).drop_duplicates()
map_df = map_df.merge(markers[["gene_upper", "gene", "group"]], on="gene_upper", how="left")
map_df.to_csv(OUT / "marker_promoter_region_map.tsv", sep="\t", index=False)

print("[INFO] marker genes found:", map_df["gene_upper"].nunique(), flush=True)
print("[INFO] marker promoter regions:", map_df["region_id"].nunique(), flush=True)

marker_regions = set(map_df["region_id"].astype(str))

print("[INFO] read P7 cluster", flush=True)
p7 = pd.read_csv(P7_UMAP, sep="\t")
p7.columns = [str(x).strip() for x in p7.columns]
cluster_col = None
for c in ["cluster", "leiden", "kmeans", "KMeans"]:
    if c in p7.columns:
        cluster_col = c
        break
if cluster_col is None:
    cand = [c for c in p7.columns if "cluster" in c.lower()]
    if cand:
        cluster_col = cand[0]
if cluster_col is None:
    raise ValueError(f"cannot find cluster column in {P7_UMAP}: {p7.columns.tolist()}")

p7["CB"] = p7["CB"].astype(str)
p7["cluster"] = pd.to_numeric(p7[cluster_col], errors="coerce")
p7 = p7[["CB", "cluster"]].dropna()
p7["cluster"] = p7["cluster"].astype(int)

print("[INFO] read selected cells", flush=True)
cell_maps = []
for version, path in CELL_FILES.items():
    c = pd.read_csv(path, sep="\t")
    c["CB"] = c["CB"].astype(str)
    c = c[["CB"]].drop_duplicates()
    c = c.merge(p7, on="CB", how="inner")
    c["version"] = version
    cell_maps.append(c)
    print(f"[INFO] {version}: {len(c)} cells", flush=True)

cell_map = pd.concat(cell_maps, ignore_index=True)
cell_set = set(cell_map["CB"])

print("[INFO] scan matrix marker rows", flush=True)
parts = []
usecols = ["CB", "region_id", "modified_calls", "covered_calls"]

for chunk in pd.read_csv(REGION_COUNTS, sep="\t", usecols=usecols, chunksize=1000000):
    chunk["CB"] = chunk["CB"].astype(str)
    chunk["region_id"] = chunk["region_id"].astype(str)

    sub = chunk[chunk["CB"].isin(cell_set) & chunk["region_id"].isin(marker_regions)].copy()
    if sub.empty:
        continue

    sub = sub.merge(map_df[["region_id", "group", "gene"]], on="region_id", how="inner")
    sub = sub.merge(cell_map, on="CB", how="inner")

    sub["modified_calls"] = pd.to_numeric(sub["modified_calls"], errors="coerce").fillna(0)
    sub["covered_calls"] = pd.to_numeric(sub["covered_calls"], errors="coerce").fillna(0)

    g = sub.groupby(["version", "CB", "cluster", "group"], as_index=False)[["modified_calls", "covered_calls"]].sum()
    parts.append(g)

if not parts:
    raise RuntimeError("no marker records found")

cell_group = pd.concat(parts, ignore_index=True)
cell_group = cell_group.groupby(["version", "CB", "cluster", "group"], as_index=False)[["modified_calls", "covered_calls"]].sum()
cell_group["modified_rate"] = cell_group["modified_calls"] / cell_group["covered_calls"].replace(0, np.nan)
cell_group["activity"] = 1 - cell_group["modified_rate"]

cell_group.to_csv(OUT / "cell_marker_group_activity.tsv", sep="\t", index=False)

# 只保留该细胞在该marker类别里有足够覆盖的记录
f = cell_group[cell_group["covered_calls"] >= MIN_GROUP_COV_PER_CELL].copy()

summary = f.groupby(["version", "cluster", "group"], as_index=False).agg(
    n_cells=("CB", "nunique"),
    median_activity=("activity", "median"),
    mean_activity=("activity", "mean"),
    median_cov=("covered_calls", "median"),
    total_cov=("covered_calls", "sum"),
    total_mod=("modified_calls", "sum"),
)
summary["weighted_activity"] = 1 - summary["total_mod"] / summary["total_cov"].replace(0, np.nan)
summary.to_csv(OUT / "cluster_marker_group_activity_summary.tsv", sep="\t", index=False)

for version, d in summary.groupby("version"):
    mat = d.pivot_table(index="group", columns="cluster", values="median_activity", aggfunc="mean")
    nmat = d.pivot_table(index="group", columns="cluster", values="n_cells", aggfunc="sum")

    mat.to_csv(OUT / f"{version}.module_median_activity.tsv", sep="\t")
    nmat.to_csv(OUT / f"{version}.module_n_cells.tsv", sep="\t")

    # 对每个marker类别做row z-score，但现在每一行是一个cell type module，不是单基因
    mat_fill = mat.apply(lambda x: x.fillna(x.mean()), axis=1)
    z = mat_fill.sub(mat_fill.mean(axis=1), axis=0)
    z = z.div(mat_fill.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    z.to_csv(OUT / f"{version}.module_median_activity_zscore.tsv", sep="\t")

    plt.figure(figsize=(7, max(4, 0.45 * z.shape[0])))
    im = plt.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    plt.xticks(range(z.shape[1]), [f"C{x}" for x in z.columns], rotation=45, ha="right")
    plt.yticks(range(z.shape[0]), z.index)
    plt.xlabel("P7 cluster")
    plt.ylabel("Marker module")
    plt.title(f"{version}: marker module activity")
    plt.colorbar(im, label="row z-score of median activity")
    plt.tight_layout()
    plt.savefig(OUT / f"{version}.module_activity_heatmap.png", dpi=240)
    plt.close()

    plt.figure(figsize=(7, max(4, 0.45 * nmat.shape[0])))
    im = plt.imshow(nmat.fillna(0).values, aspect="auto")
    plt.xticks(range(nmat.shape[1]), [f"C{x}" for x in nmat.columns], rotation=45, ha="right")
    plt.yticks(range(nmat.shape[0]), nmat.index)
    plt.xlabel("P7 cluster")
    plt.ylabel("Marker module")
    plt.title(f"{version}: cells with marker-module coverage")
    plt.colorbar(im, label="n_cells with group coverage")
    plt.tight_layout()
    plt.savefig(OUT / f"{version}.module_coverage_ncells_heatmap.png", dpi=240)
    plt.close()

print("[INFO] done", flush=True)
