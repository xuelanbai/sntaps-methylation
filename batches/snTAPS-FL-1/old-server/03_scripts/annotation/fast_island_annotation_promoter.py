import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.cluster import DBSCAN
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy.stats import fisher_exact
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False

base = Path("/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation")

UMAP = base / "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"

MATRIX = base / "08_qc_summary/final_single_sample/region_matrix_covfilter/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.cov10.region_counts.tsv.gz"

OUT = base / "08_qc_summary/final_single_sample/annotation_fast/P7_promoter_cov10_island_marker"
OUT.mkdir(parents=True, exist_ok=True)

# DBSCAN参数：先把明显分离的小岛分开
EPS = 0.45
MIN_SAMPLES = 40
MIN_CLUSTER_CELLS = 500
MIN_OBS_IN = 100
MIN_OBS_OUT = 300
MIN_COV_IN = 1000
MIN_COV_OUT = 3000
TOP_N = 100

print("[INFO] read UMAP:", UMAP)
df = pd.read_csv(UMAP, sep="\t")

if "CB" not in df.columns:
    df = df.rename(columns={df.columns[0]: "CB"})

need_cols = ["CB", "UMAP1", "UMAP2"]
for c in need_cols:
    if c not in df.columns:
        raise ValueError(f"missing column: {c}")

xy = df[["UMAP1", "UMAP2"]].values
raw_lab = DBSCAN(eps=EPS, min_samples=MIN_SAMPLES, n_jobs=8).fit_predict(xy)

# 按cluster大小重新编号，最大岛为0；noise保留为noise
vc_raw = pd.Series(raw_lab).value_counts()
valid_labs = [x for x in vc_raw.index if x != -1]
valid_labs = sorted(valid_labs, key=lambda x: vc_raw[x], reverse=True)
lab_map = {old: str(i) for i, old in enumerate(valid_labs)}

df["island"] = [lab_map.get(x, "noise") for x in raw_lab]
df.to_csv(OUT / "P7_promoter_cov10.DBSCAN_island_umap.tsv", sep="\t", index=False)

counts = df["island"].value_counts().rename_axis("island").reset_index(name="n_cells")
counts["fraction"] = counts["n_cells"] / counts["n_cells"].sum()
counts.to_csv(OUT / "P7_promoter_cov10.DBSCAN_island_counts.tsv", sep="\t", index=False)

# QC汇总
qc_cols = [c for c in ["reads_used", "covered_calls", "covered_bins", "modified_rate"] if c in df.columns]
if qc_cols:
    qc = df[df["island"] != "noise"].groupby("island")[qc_cols].agg(["count", "median", "mean"])
    qc.to_csv(OUT / "P7_promoter_cov10.DBSCAN_island_qc_summary.tsv", sep="\t")

# 画岛级分群图
plot_df = df.copy()
fig, ax = plt.subplots(figsize=(5.5, 5.0))
labs = [x for x in counts["island"].tolist() if x != "noise"]
cmap = plt.get_cmap("tab20")

# 先画noise
sub = plot_df[plot_df["island"] == "noise"]
if len(sub) > 0:
    ax.scatter(sub["UMAP1"], sub["UMAP2"], s=0.4, color="lightgrey", alpha=0.35, linewidths=0, label="noise")

for i, lab in enumerate(labs):
    sub = plot_df[plot_df["island"] == lab]
    ax.scatter(sub["UMAP1"], sub["UMAP2"], s=0.55, color=cmap(i % 20), alpha=0.78, linewidths=0, label=lab)

ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.set_title("Promoter residual UMAP")
if len(labs) <= 25:
    ax.legend(markerscale=5, fontsize=5, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
fig.tight_layout()
fig.savefig(OUT / "P7_promoter_cov10.DBSCAN_island_umap_by_island.png", dpi=600)
plt.close(fig)

# density图
fig, ax = plt.subplots(figsize=(5.5, 5.0))
hb = ax.hexbin(df["UMAP1"], df["UMAP2"], gridsize=90, mincnt=1)
cb = plt.colorbar(hb, ax=ax, fraction=0.046, pad=0.04)
cb.set_label("Cell density")
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.set_title("Promoter residual UMAP density")
fig.tight_layout()
fig.savefig(OUT / "P7_promoter_cov10.DBSCAN_island_density.png", dpi=600)
plt.close(fig)

# 只给足够大的岛做marker
use_islands = counts[(counts["island"] != "noise") & (counts["n_cells"] >= MIN_CLUSTER_CELLS)]["island"].astype(str).tolist()
print("[INFO] islands for marker:", use_islands)

label_map2 = df[df["island"].isin(use_islands)].set_index("CB")["island"].to_dict()
if not label_map2:
    raise ValueError("No valid islands for marker analysis.")

tmp = OUT / "partial_cluster_region_agg.tsv"
if tmp.exists():
    tmp.unlink()

print("[INFO] scan matrix:", MATRIX)
reader = pd.read_csv(MATRIX, sep="\t", chunksize=1_500_000)

first = True
region_col = None

for i, chunk in enumerate(reader, 1):
    if region_col is None:
        if "region_id" in chunk.columns:
            region_col = "region_id"
        elif "bin_id" in chunk.columns:
            region_col = "bin_id"
        else:
            raise ValueError("Cannot find region_id/bin_id column in matrix.")

    if "CB" not in chunk.columns:
        raise ValueError("Matrix missing CB column.")

    chunk = chunk[chunk["CB"].isin(label_map2)]
    if chunk.empty:
        continue

    chunk["island"] = chunk["CB"].map(label_map2)
    g = chunk.groupby(["island", region_col], observed=True).agg(
        modified_calls=("modified_calls", "sum"),
        covered_calls=("covered_calls", "sum"),
        n_obs=("CB", "size")
    ).reset_index()

    g.to_csv(tmp, sep="\t", index=False, mode="a", header=first)
    first = False

    if i % 10 == 0:
        print("[INFO] processed chunks:", i)

print("[INFO] aggregate partial table")
part = pd.read_csv(tmp, sep="\t")
agg = part.groupby(["island", region_col], observed=True).agg(
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

# 低甲基化promoter：更可能对应该cluster活跃marker
hypo = (
    m.sort_values(["island", "delta_rate"], ascending=[True, True])
     .groupby("island", observed=True)
     .head(TOP_N)
     .copy()
)

# 高甲基化promoter：辅助参考
hyper = (
    m.sort_values(["island", "delta_rate"], ascending=[True, False])
     .groupby("island", observed=True)
     .head(TOP_N)
     .copy()
)

def add_pval(tab):
    if not HAS_SCIPY:
        tab["fisher_pvalue"] = np.nan
        return tab
    pvals = []
    for _, r in tab.iterrows():
        a = int(r["modified_calls"])
        b = int(r["covered_calls"] - r["modified_calls"])
        c = int(r["out_modified"])
        d = int(r["out_covered"] - r["out_modified"])
        try:
            _, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        except Exception:
            p = np.nan
        pvals.append(p)
    tab["fisher_pvalue"] = pvals
    return tab

hypo = add_pval(hypo)
hyper = add_pval(hyper)

hypo.to_csv(OUT / "P7_promoter_cov10.island_hypomethylated_promoter_markers.top100.tsv", sep="\t", index=False)
hyper.to_csv(OUT / "P7_promoter_cov10.island_hypermethylated_promoter_markers.top100.tsv", sep="\t", index=False)

# 常见细胞类型marker，用低甲基化promoter初步匹配
marker_panel = {
    "T_cell": ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "IL7R", "CD2", "CD247"],
    "NK_cell": ["NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "KLRK1"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD74", "BANK1", "CD19"],
    "Plasma_cell": ["MZB1", "JCHAIN", "XBP1", "IGHG1", "IGKC"],
    "Monocyte_macrophage": ["LYZ", "LST1", "TYROBP", "FCGR3A", "CD14", "CSF1R", "C1QA", "C1QB", "C1QC"],
    "Neutrophil_granulocyte": ["S100A8", "S100A9", "MPO", "ELANE", "CSF3R"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "FLT1", "EMCN", "CDH5", "CLDN5"],
    "Fibroblast_stromal": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"],
    "Smooth_muscle_pericyte": ["ACTA2", "MYH11", "TAGLN", "RGS5", "PDGFRB", "CSPG4"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7"],
    "Muscle": ["MYOD1", "MYOG", "MYH1", "MYH2", "MYH3", "ACTA1", "TNNT1", "TNNT3", "DES", "MYL1"],
    "Hepatocyte": ["ALB", "APOA1", "APOB", "TTR", "CYP3A4", "CYP2E1", "ASGR1"],
    "Erythroid": ["HBA1", "HBA2", "HBB", "ALAS2", "GYPA"]
}

hypo["region_text"] = hypo[region_col].astype(str).str.upper()

ann_rows = []
hit_rows = []

for island, sub in hypo.groupby("island", observed=True):
    sub = sub.copy()
    sub["rank"] = np.arange(1, len(sub) + 1)

    for celltype, genes in marker_panel.items():
        hits = []
        score = 0.0
        for gene in genes:
            mask = sub["region_text"].str.contains(gene.upper(), regex=False, na=False)
            if mask.any():
                best = sub.loc[mask].sort_values("rank").iloc[0]
                hits.append(f"{gene}:{best[region_col]}:rank{int(best['rank'])}:delta{best['delta_rate']:.4f}")
                score += 1.0 / np.sqrt(float(best["rank"]))

        if hits:
            hit_rows.append({
                "island": island,
                "candidate_celltype": celltype,
                "score": score,
                "n_marker_hits": len(hits),
                "hits": ";".join(hits)
            })

if hit_rows:
    hits_df = pd.DataFrame(hit_rows).sort_values(["island", "score"], ascending=[True, False])
    hits_df.to_csv(OUT / "P7_promoter_cov10.island_marker_panel_hits.tsv", sep="\t", index=False)

    best = hits_df.groupby("island", observed=True).head(3)
    best.to_csv(OUT / "P7_promoter_cov10.island_annotation_suggestion.top3.tsv", sep="\t", index=False)
else:
    pd.DataFrame(columns=["island", "candidate_celltype", "score", "n_marker_hits", "hits"]).to_csv(
        OUT / "P7_promoter_cov10.island_marker_panel_hits.tsv", sep="\t", index=False
    )

print("[DONE]", OUT)
print("[MAIN OUTPUTS]")
for f in [
    "P7_promoter_cov10.DBSCAN_island_counts.tsv",
    "P7_promoter_cov10.DBSCAN_island_umap_by_island.png",
    "P7_promoter_cov10.DBSCAN_island_density.png",
    "P7_promoter_cov10.island_hypomethylated_promoter_markers.top100.tsv",
    "P7_promoter_cov10.island_marker_panel_hits.tsv",
    "P7_promoter_cov10.island_annotation_suggestion.top3.tsv",
]:
    print(OUT / f)
