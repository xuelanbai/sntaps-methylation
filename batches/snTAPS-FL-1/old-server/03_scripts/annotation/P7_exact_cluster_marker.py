import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import fisher_exact

BASE = Path("/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation")

UMAP = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"

MATRIX = BASE / "08_qc_summary/final_single_sample/region_matrix_covfilter/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.cov10.region_counts.tsv.gz"

OUT = BASE / "08_qc_summary/final_single_sample/annotation_fast/P7_exact_cluster_marker"
OUT.mkdir(parents=True, exist_ok=True)

MIN_CLUSTER_CELLS = 300
MIN_OBS_IN = 100
MIN_OBS_OUT = 300
MIN_COV_IN = 1000
MIN_COV_OUT = 3000
TOP_N = 300

marker_panel = {
    "T_cell": ["CD3D","CD3E","CD3G","TRAC","TRBC1","TRBC2","IL7R","CD2","CD247","LCK"],
    "NK_cell": ["NKG7","GNLY","PRF1","GZMB","KLRD1","KLRK1"],
    "B_cell": ["MS4A1","CD79A","CD79B","CD74","BANK1","CD19","PAX5"],
    "Plasma_cell": ["MZB1","JCHAIN","XBP1","IGHG1","IGKC"],
    "Myeloid": ["LYZ","LST1","TYROBP","CD14","CD68","CSF1R","C1QA","C1QB","C1QC","FCGR3A"],
    "Neutrophil": ["S100A8","S100A9","MPO","ELANE","CSF3R"],
    "Endothelial": ["PECAM1","VWF","KDR","FLT1","EMCN","CDH5","CLDN5"],
    "Fibroblast_stromal": ["COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA"],
    "Smooth_muscle_pericyte": ["ACTA2","MYH11","TAGLN","RGS5","PDGFRB","CSPG4"],
    "Epithelial": ["EPCAM","KRT8","KRT18","KRT19","KRT7"],
    "Muscle": ["MYOD1","MYOG","MYH1","MYH2","MYH3","ACTA1","TNNT1","TNNT3","DES","MYL1","MYL3"],
    "Hepatocyte": ["ALB","AFP","APOA1","APOB","TTR","ASGR1","CYP3A4","CYP2E1","SERPINA1"],
    "Erythroid": ["HBA1","HBA2","HBB","ALAS2","GYPA","KLF1","GATA1","EPOR"],
}

def get_symbol(region_id):
    s = str(region_id)
    if "|" in s:
        return s.split("|")[-1].upper()
    return ""

print("[INFO] read UMAP")
meta = pd.read_csv(UMAP, sep="\t")

if "CB" not in meta.columns:
    meta = meta.rename(columns={meta.columns[0]: "CB"})

cluster_col = None
for c in ["cluster", "leiden", "kmeans", "Cluster"]:
    if c in meta.columns:
        cluster_col = c
        break

if cluster_col is None:
    raise ValueError("UMAP文件里没有cluster/leiden/kmeans列，请先确认P7的umap.tsv是否带cluster。")

meta["cluster"] = meta[cluster_col].astype(str)

cluster_counts = meta["cluster"].value_counts().rename_axis("cluster").reset_index(name="n_cells")
cluster_counts["fraction"] = cluster_counts["n_cells"] / cluster_counts["n_cells"].sum()
cluster_counts.to_csv(OUT / "P7.cluster_counts.tsv", sep="\t", index=False)

use_clusters = cluster_counts.loc[cluster_counts["n_cells"] >= MIN_CLUSTER_CELLS, "cluster"].astype(str).tolist()
label_map = meta[meta["cluster"].isin(use_clusters)].set_index("CB")["cluster"].to_dict()

print("[INFO] clusters used:", use_clusters)
print("[INFO] cells used:", len(label_map))

tmp = OUT / "P7.partial_cluster_region_agg.tsv"
if tmp.exists():
    tmp.unlink()

reader = pd.read_csv(MATRIX, sep="\t", chunksize=1500000)
first = True
region_col = None

for i, chunk in enumerate(reader, 1):
    if region_col is None:
        if "region_id" in chunk.columns:
            region_col = "region_id"
        elif "bin_id" in chunk.columns:
            region_col = "bin_id"
        else:
            raise ValueError("找不到region_id/bin_id列")

    chunk = chunk[chunk["CB"].isin(label_map)]
    if chunk.empty:
        continue

    chunk["cluster"] = chunk["CB"].map(label_map)

    g = chunk.groupby(["cluster", region_col], observed=True).agg(
        modified_calls=("modified_calls", "sum"),
        covered_calls=("covered_calls", "sum"),
        n_obs=("CB", "size")
    ).reset_index()

    g.to_csv(tmp, sep="\t", index=False, mode="a", header=first)
    first = False

    if i % 10 == 0:
        print("[INFO] chunks processed:", i)

print("[INFO] aggregate")
part = pd.read_csv(tmp, sep="\t")

agg = part.groupby(["cluster", region_col], observed=True).agg(
    modified_calls=("modified_calls", "sum"),
    covered_calls=("covered_calls", "sum"),
    n_obs=("n_obs", "sum")
).reset_index()

total = agg.groupby(region_col, observed=True).agg(
    total_modified=("modified_calls", "sum"),
    total_covered=("covered_calls", "sum"),
    total_obs=("n_obs", "sum")
).reset_index()

m = agg.merge(total, on=region_col, how="left")
m["out_modified"] = m["total_modified"] - m["modified_calls"]
m["out_covered"] = m["total_covered"] - m["covered_calls"]
m["out_obs"] = m["total_obs"] - m["n_obs"]

m = m[
    (m["covered_calls"] >= MIN_COV_IN) &
    (m["out_covered"] >= MIN_COV_OUT) &
    (m["n_obs"] >= MIN_OBS_IN) &
    (m["out_obs"] >= MIN_OBS_OUT)
].copy()

m["rate_in"] = m["modified_calls"] / m["covered_calls"]
m["rate_out"] = m["out_modified"] / m["out_covered"]
m["delta_rate"] = m["rate_in"] - m["rate_out"]
m["abs_delta"] = m["delta_rate"].abs()
m["gene_symbol"] = m[region_col].map(get_symbol)

pvals = []
for _, r in m.iterrows():
    a = int(r["modified_calls"])
    b = int(r["covered_calls"] - r["modified_calls"])
    c = int(r["out_modified"])
    d = int(r["out_covered"] - r["out_modified"])
    try:
        _, p = fisher_exact([[a,b],[c,d]], alternative="two-sided")
    except Exception:
        p = np.nan
    pvals.append(p)

m["fisher_pvalue"] = pvals

# 低甲基化promoter作为候选活跃marker
hypo = (
    m.sort_values(["cluster", "delta_rate"], ascending=[True, True])
     .groupby("cluster", observed=True)
     .head(TOP_N)
     .copy()
)

hypo.to_csv(OUT / "P7.cluster_hypomethylated_promoter_markers.top300.tsv", sep="\t", index=False)

# 高甲基化也输出，备用
hyper = (
    m.sort_values(["cluster", "delta_rate"], ascending=[True, False])
     .groupby("cluster", observed=True)
     .head(TOP_N)
     .copy()
)
hyper.to_csv(OUT / "P7.cluster_hypermethylated_promoter_markers.top300.tsv", sep="\t", index=False)

# 精确marker匹配，不再用contains
hit_rows = []

for cluster, sub in hypo.groupby("cluster", observed=True):
    sub = sub.copy()
    sub["rank"] = np.arange(1, len(sub) + 1)
    sym_to_rows = {g: x for g, x in sub.groupby("gene_symbol") if g != ""}

    for celltype, genes in marker_panel.items():
        hits = []
        score = 0.0
        for gene in genes:
            gene_u = gene.upper()
            if gene_u in sym_to_rows:
                best = sym_to_rows[gene_u].sort_values("rank").iloc[0]
                rank = int(best["rank"])
                delta = float(best["delta_rate"])
                hits.append(f"{gene}:{best[region_col]}:rank{rank}:delta{delta:.4f}")
                score += 1 / np.sqrt(rank)

        if hits:
            hit_rows.append({
                "cluster": cluster,
                "candidate_celltype": celltype,
                "score": score,
                "n_marker_hits": len(hits),
                "hits": ";".join(hits)
            })

hits_df = pd.DataFrame(hit_rows)

if len(hits_df) > 0:
    hits_df = hits_df.sort_values(["cluster", "score"], ascending=[True, False])
    hits_df.to_csv(OUT / "P7.cluster_marker_panel_hits.EXACT.tsv", sep="\t", index=False)

    suggestion = hits_df.groupby("cluster", observed=True).head(5)
    suggestion.to_csv(OUT / "P7.cluster_annotation_suggestion.EXACT.top5.tsv", sep="\t", index=False)
else:
    pd.DataFrame(columns=["cluster","candidate_celltype","score","n_marker_hits","hits"]).to_csv(
        OUT / "P7.cluster_marker_panel_hits.EXACT.tsv", sep="\t", index=False
    )
    pd.DataFrame(columns=["cluster","candidate_celltype","score","n_marker_hits","hits"]).to_csv(
        OUT / "P7.cluster_annotation_suggestion.EXACT.top5.tsv", sep="\t", index=False
    )

print("[DONE]", OUT)
