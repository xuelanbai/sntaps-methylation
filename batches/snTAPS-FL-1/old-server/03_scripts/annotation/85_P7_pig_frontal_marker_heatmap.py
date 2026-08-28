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

# 诊断已经证明cov10有64098条marker x selected-cell记录，优先用cov10
REGION_COUNTS = BASE / "08_qc_summary/final_single_sample/region_matrix_covfilter/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.cov10.region_counts.tsv.gz"

CELL_FILES = {
    "P7_QCpass_strict": BASE / "08_qc_summary/final_single_sample/cell_filter_P7_QCpass/snTAPS-FL-1.P7_QCpass_strict.global_per_cb.tsv",
    "P7_QCpass_no01_reference": BASE / "08_qc_summary/final_single_sample/cell_filter_P7_QCpass/snTAPS-FL-1.P7_QCpass_no01_reference.global_per_cb.tsv",
    "P7_QCpass_no016_reference": BASE / "08_qc_summary/final_single_sample/cell_filter_P7_QCpass/snTAPS-FL-1.P7_QCpass_no016_reference.global_per_cb.tsv",
}

P7_UMAP = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"

OUT = BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/P7_QCpass_marker"
OUT.mkdir(parents=True, exist_ok=True)

print("[INFO] read marker list", flush=True)
markers = pd.read_csv(MARKER_TSV, sep="\t")
markers.columns = [str(x).strip() for x in markers.columns]
if "gene" not in markers.columns or "group" not in markers.columns:
    raise ValueError(f"marker file columns are {markers.columns.tolist()}, expected group and gene")

markers["gene"] = markers["gene"].astype(str).str.strip()
markers["group"] = markers["group"].astype(str).str.strip()
markers = markers[(markers["gene"] != "") & (markers["gene"].str.lower() != "nan")].copy()
markers["gene_upper"] = markers["gene"].str.upper()
marker_set = set(markers["gene_upper"])
print(f"[INFO] marker genes: {len(marker_set)}", flush=True)

print("[INFO] parse promoter BED", flush=True)
bed_rows = []
with open(BED) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = re.split(r"\s+", line.strip())

        ens = [x for x in parts if re.match(r"^ENSSSCG\d+", x)]
        if not ens:
            continue
        region_id = ens[0]

        hits = []
        for x in parts:
            xu = x.upper()
            if xu in marker_set:
                hits.append(xu)

        for gene_upper in sorted(set(hits)):
            bed_rows.append((region_id, gene_upper))

map_df = pd.DataFrame(bed_rows, columns=["region_id", "gene_upper"]).drop_duplicates()
map_df = map_df.merge(markers, on="gene_upper", how="left")
map_df = map_df[["region_id", "gene", "gene_upper", "group"]].drop_duplicates()

map_df.to_csv(OUT / "marker_promoter_region_map.tsv", sep="\t", index=False)
pd.DataFrame({"found_marker_gene": sorted(map_df["gene_upper"].unique())}).to_csv(
    OUT / "found_marker_genes.tsv", sep="\t", index=False
)
pd.DataFrame({"missing_marker_gene": sorted(marker_set - set(map_df["gene_upper"]))}).to_csv(
    OUT / "missing_marker_genes.tsv", sep="\t", index=False
)

print(f"[INFO] marker regions from BED: {map_df.shape}", flush=True)
print(f"[INFO] found marker genes: {map_df['gene_upper'].nunique()}", flush=True)

marker_regions = set(map_df["region_id"].astype(str))
if len(marker_regions) == 0:
    raise RuntimeError("No marker promoter regions found from BED.")

print("[INFO] read P7 cluster map", flush=True)
p7 = pd.read_csv(P7_UMAP, sep="\t")
p7.columns = [str(x).strip() for x in p7.columns]

cluster_col = None
for c in ["cluster", "leiden", "kmeans", "KMeans", "old_cluster"]:
    if c in p7.columns:
        cluster_col = c
        break
if cluster_col is None:
    # 自动找含cluster的列
    cand = [c for c in p7.columns if "cluster" in c.lower()]
    if cand:
        cluster_col = cand[0]
if cluster_col is None:
    raise ValueError(f"cannot find cluster column in {P7_UMAP}; columns={p7.columns.tolist()}")

p7["CB"] = p7["CB"].astype(str)
p7["cluster"] = pd.to_numeric(p7[cluster_col], errors="coerce")
p7 = p7[["CB", "cluster"]].dropna().copy()
p7["cluster"] = p7["cluster"].astype(int)
print(f"[INFO] P7 cluster cells: {len(p7)}", flush=True)

print("[INFO] read selected cell lists", flush=True)
cell_maps = []
for version, path in CELL_FILES.items():
    c = pd.read_csv(path, sep="\t")
    c.columns = [str(x).strip() for x in c.columns]
    if "CB" not in c.columns:
        raise ValueError(f"CB not found in {path}; columns={c.columns.tolist()}")

    c["CB"] = c["CB"].astype(str)
    c = c[["CB"]].drop_duplicates()
    c = c.merge(p7, on="CB", how="inner")
    c["version"] = version
    cell_maps.append(c)
    print(f"[INFO] {version}: selected={len(c)}", flush=True)

cell_map = pd.concat(cell_maps, ignore_index=True)
cell_map.to_csv(OUT / "all_versions_cell_cluster_map.tsv", sep="\t", index=False)

cell_set = set(cell_map["CB"])
print(f"[INFO] total version-cell records: {len(cell_map)}", flush=True)
print(f"[INFO] unique selected CBs: {len(cell_set)}", flush=True)

print("[INFO] scan promoter matrix", flush=True)
agg_parts = []
usecols = ["CB", "region_id", "modified_calls", "covered_calls"]

total_marker_rows = 0
total_marker_cell_rows = 0

for chunk in pd.read_csv(REGION_COUNTS, sep="\t", usecols=usecols, chunksize=1000000):
    chunk["CB"] = chunk["CB"].astype(str)
    chunk["region_id"] = chunk["region_id"].astype(str)

    mr = chunk["region_id"].isin(marker_regions)
    total_marker_rows += int(mr.sum())

    sub = chunk[mr & chunk["CB"].isin(cell_set)].copy()
    total_marker_cell_rows += len(sub)

    if sub.empty:
        continue

    sub["modified_calls"] = pd.to_numeric(sub["modified_calls"], errors="coerce").fillna(0)
    sub["covered_calls"] = pd.to_numeric(sub["covered_calls"], errors="coerce").fillna(0)

    sub = sub.merge(cell_map, on="CB", how="inner")
    sub = sub.merge(map_df, on="region_id", how="inner")

    g = sub.groupby(["version", "cluster", "group", "gene"], as_index=False)[["modified_calls", "covered_calls"]].sum()
    agg_parts.append(g)

print(f"[INFO] total rows with marker region: {total_marker_rows}", flush=True)
print(f"[INFO] total rows with marker region AND selected cells: {total_marker_cell_rows}", flush=True)

if not agg_parts:
    raise RuntimeError("No marker records found after merge. Check cell_map and region_id.")

agg = pd.concat(agg_parts, ignore_index=True)
agg = agg.groupby(["version", "cluster", "group", "gene"], as_index=False)[["modified_calls", "covered_calls"]].sum()

agg["modified_rate"] = agg["modified_calls"] / agg["covered_calls"].replace(0, np.nan)
agg["promoter_activity"] = 1 - agg["modified_rate"]
agg.to_csv(OUT / "P7_marker_cluster_activity_long.tsv", sep="\t", index=False)

print("[INFO] plot heatmaps", flush=True)
for version, d in agg.groupby("version"):
    d = d[d["covered_calls"] > 0].copy()
    d["label"] = d["group"] + ":" + d["gene"]

    order = markers.copy()
    order["label"] = order["group"] + ":" + order["gene"]

    mat = d.pivot_table(index="label", columns="cluster", values="promoter_activity", aggfunc="mean")
    keep_order = [x for x in order["label"].tolist() if x in mat.index]
    mat = mat.loc[keep_order]
    mat = mat.dropna(how="all")

    mat_fill = mat.apply(lambda x: x.fillna(x.mean()), axis=1)
    z = mat_fill.sub(mat_fill.mean(axis=1), axis=0)
    sd = mat_fill.std(axis=1).replace(0, np.nan)
    z = z.div(sd, axis=0).fillna(0)

    mat.to_csv(OUT / f"{version}.marker_promoter_activity.tsv", sep="\t")
    z.to_csv(OUT / f"{version}.marker_promoter_activity_zscore.tsv", sep="\t")

    top_rows = []
    for cl in z.columns:
        tmp = z[[cl]].sort_values(cl, ascending=False).head(20).copy()
        tmp["cluster"] = cl
        tmp["marker"] = tmp.index
        tmp["zscore"] = tmp[cl].values
        top_rows.append(tmp[["cluster", "marker", "zscore"]])
    pd.concat(top_rows).to_csv(OUT / f"{version}.top_marker_by_cluster.tsv", sep="\t", index=False)

    h = max(7, 0.25 * z.shape[0])
    w = max(7, 0.6 * z.shape[1] + 3)
    plt.figure(figsize=(w, h))
    im = plt.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    plt.xticks(range(z.shape[1]), [f"C{x}" for x in z.columns], rotation=45, ha="right")
    plt.yticks(range(z.shape[0]), z.index, fontsize=7)
    plt.xlabel("P7 cluster")
    plt.ylabel("Pig frontal lobe marker")
    plt.title(f"{version}: promoter hypomethylation marker activity")
    plt.colorbar(im, label="row z-score of promoter activity")
    plt.tight_layout()
    plt.savefig(OUT / f"{version}.pig_frontal_marker_heatmap.png", dpi=240)
    plt.close()

print("[INFO] done", flush=True)
