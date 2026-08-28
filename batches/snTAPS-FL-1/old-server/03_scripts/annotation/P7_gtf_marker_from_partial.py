import re
import gzip
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import fisher_exact
import sys

gtf = Path(sys.argv[1])

base = Path("/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation")
partial = base / "08_qc_summary/final_single_sample/annotation_fast/P7_exact_cluster_marker/P7.partial_cluster_region_agg.tsv"
out = base / "08_qc_summary/final_single_sample/annotation_fast/P7_gtf_marker"
out.mkdir(parents=True, exist_ok=True)

TOP_N = 1000
MIN_OBS_IN = 100
MIN_OBS_OUT = 300
MIN_COV_IN = 1000
MIN_COV_OUT = 3000

marker_panel = {
    "T_cell": ["CD3D","CD3E","CD3G","TRAC","TRBC1","TRBC2","IL7R","CD2","CD247","LCK"],
    "NK_cell": ["NKG7","GNLY","PRF1","GZMB","KLRD1","KLRK1"],
    "B_cell": ["MS4A1","CD79A","CD79B","CD74","BANK1","CD19","PAX5"],
    "Myeloid": ["LYZ","LST1","TYROBP","CD14","CD68","CSF1R","C1QA","C1QB","C1QC","FCGR3A"],
    "Neutrophil": ["S100A8","S100A9","MPO","ELANE","CSF3R"],
    "Endothelial": ["PECAM1","VWF","KDR","FLT1","EMCN","CDH5","CLDN5"],
    "Fibroblast_stromal": ["COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA"],
    "Smooth_muscle_pericyte": ["ACTA2","MYH11","TAGLN","RGS5","PDGFRB","CSPG4"],
    "Skeletal_muscle": ["MYOD1","MYOG","MYF5","MYF6","PAX7","ACTA1","DES","TNNT1","TNNT3","MYL1","MYL3","MYH1","MYH2","MYH3"],
    "Adipocyte": ["ADIPOQ","PPARG","FABP4","PLIN1","CEBPA","LPL"],
    "Epithelial": ["EPCAM","KRT8","KRT18","KRT19","KRT7"],
    "Hepatocyte": ["ALB","AFP","APOA1","APOB","TTR","ASGR1","CYP3A4","CYP2E1","SERPINA1"],
    "Erythroid": ["HBA1","HBA2","HBB","ALAS2","GYPA","KLF1","GATA1"],
}

def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "r")

def parse_attrs(attr):
    d = {}
    # GTF: gene_id "xxx"; gene_name "yyy";
    for k, v in re.findall(r'(\S+)\s+"([^"]+)"', attr):
        d[k] = v
    # GFF3: ID=xxx;Name=yyy
    for item in attr.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

gene_map = {}
with opener(gtf) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        feat = parts[2]
        if feat not in ["gene", "transcript", "mRNA"]:
            continue
        attrs = parse_attrs(parts[8])
        gid = attrs.get("gene_id") or attrs.get("ID") or attrs.get("Parent")
        name = attrs.get("gene_name") or attrs.get("gene_symbol") or attrs.get("Name") or attrs.get("gene") or gid
        if gid:
            gid0 = gid.split(".")[0]
            gene_map[gid0] = name.upper()

pd.DataFrame({"gene_id": list(gene_map.keys()), "gene_symbol": list(gene_map.values())}).to_csv(
    out / "pig_gene_id_to_symbol.tsv", sep="\t", index=False
)

print("[INFO] gene map size:", len(gene_map))

part = pd.read_csv(partial, sep="\t")
region_col = "region_id" if "region_id" in part.columns else "bin_id"

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

m["gene_id"] = m[region_col].astype(str).str.split("|").str[0].str.split(".").str[0]
m["gene_symbol"] = m[region_col].astype(str).str.split("|").str[-1]
m.loc[~m["gene_symbol"].str.match(r"^[A-Za-z]", na=False), "gene_symbol"] = ""
m["mapped_symbol"] = m["gene_symbol"].replace("", np.nan)
m["mapped_symbol"] = m["mapped_symbol"].fillna(m["gene_id"].map(gene_map)).fillna("")
m["mapped_symbol"] = m["mapped_symbol"].str.upper()

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

hypo = (
    m.sort_values(["cluster", "delta_rate"], ascending=[True, True])
     .groupby("cluster", observed=True)
     .head(TOP_N)
     .copy()
)

hypo.to_csv(out / "P7.cluster_hypomethylated_promoter_markers.with_gtf.top1000.tsv", sep="\t", index=False)

hit_rows = []
for cluster, sub in hypo.groupby("cluster", observed=True):
    sub = sub.copy()
    sub["rank"] = np.arange(1, len(sub) + 1)

    for celltype, genes in marker_panel.items():
        hits = []
        score = 0.0
        for gene in genes:
            gene_u = gene.upper()
            z = sub[sub["mapped_symbol"] == gene_u]
            if len(z) > 0:
                best = z.iloc[0]
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

hits = pd.DataFrame(hit_rows)
if len(hits) > 0:
    hits = hits.sort_values(["cluster", "score", "n_marker_hits"], ascending=[True, False, False])
    hits.to_csv(out / "P7.cluster_annotation_suggestion.with_gtf.top5.tsv", sep="\t", index=False)
else:
    hits = pd.DataFrame(columns=["cluster","candidate_celltype","score","n_marker_hits","hits"])
    hits.to_csv(out / "P7.cluster_annotation_suggestion.with_gtf.top5.tsv", sep="\t", index=False)

# 每个cluster输出前30个有symbol的marker，方便人工判断
sym_top = hypo[hypo["mapped_symbol"] != ""].copy()
sym_top = sym_top.sort_values(["cluster", "delta_rate"], ascending=[True, True]).groupby("cluster", observed=True).head(30)
sym_top.to_csv(out / "P7.cluster_top30_symbol_markers.tsv", sep="\t", index=False)

print("[DONE]", out)
