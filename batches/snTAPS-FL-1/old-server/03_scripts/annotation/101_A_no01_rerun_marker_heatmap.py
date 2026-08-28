#!/usr/bin/env python3
from pathlib import Path
import os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A_UMAP = Path(os.environ.get(
    "A_UMAP",
    "08_qc_summary/final_single_sample/region_umap_paper/P7_no01_promoter_residual_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"
))

BED = Path("02_reference/pig/gene_regions/pig_promoter_TSS_1.5kb.bed")

MATRIX1 = Path("08_qc_summary/final_single_sample/region_matrix_P7_no01/snTAPS-FL-1.P7_no01.promoter_TSS_1.5kb.region_counts.tsv.gz")
MATRIX2 = Path("08_qc_summary/final_single_sample/region_matrix/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.region_counts.tsv.gz")
MATRIX = MATRIX2

OUTDIR = Path("08_qc_summary/final_single_sample/marker_pig_frontal_lobe/A_no01_rerun_marker")
OUTDIR.mkdir(parents=True, exist_ok=True)

markers = [
    ("Excitatory_neuron","RBFOX3"),("Excitatory_neuron","SNAP25"),("Excitatory_neuron","SYT1"),
    ("Excitatory_neuron","MAP2"),("Excitatory_neuron","SLC17A7"),("Excitatory_neuron","CAMK2A"),
    ("Excitatory_neuron","NEFL"),
    ("Upper_layer_excitatory","SATB2"),("Upper_layer_excitatory","CUX1"),
    ("Upper_layer_excitatory","CUX2"),("Upper_layer_excitatory","RORB"),
    ("Deep_layer_excitatory","BCL11B"),("Deep_layer_excitatory","FEZF2"),
    ("Deep_layer_excitatory","TBR1"),("Deep_layer_excitatory","FOXP2"),
    ("Inhibitory_neuron","GAD1"),("Inhibitory_neuron","GAD2"),("Inhibitory_neuron","SLC32A1"),
    ("Inhibitory_neuron","PVALB"),("Inhibitory_neuron","SST"),("Inhibitory_neuron","VIP"),
    ("Inhibitory_neuron","LAMP5"),("Inhibitory_neuron","RELN"),("Inhibitory_neuron","NPY"),
    ("Astrocyte","AQP4"),("Astrocyte","GFAP"),("Astrocyte","ALDH1L1"),
    ("Astrocyte","SLC1A2"),("Astrocyte","SLC1A3"),("Astrocyte","GJA1"),("Astrocyte","SOX9"),
    ("Oligodendrocyte","MBP"),("Oligodendrocyte","MOG"),("Oligodendrocyte","PLP1"),
    ("Oligodendrocyte","MAG"),("Oligodendrocyte","CLDN11"),
    ("OPC","PDGFRA"),("OPC","CSPG4"),("OPC","VCAN"),("OPC","OLIG1"),("OPC","OLIG2"),("OPC","SOX10"),
    ("Microglia","AIF1"),("Microglia","CX3CR1"),("Microglia","P2RY12"),
    ("Microglia","CSF1R"),("Microglia","TYROBP"),("Microglia","TREM2"),
    ("Microglia","C1QA"),("Microglia","C1QB"),("Microglia","C1QC"),
    ("Endothelial","PECAM1"),("Endothelial","CLDN5"),("Endothelial","VWF"),("Endothelial","KDR"),("Endothelial","FLT1"),
    ("Pericyte_SMC","PDGFRB"),("Pericyte_SMC","RGS5"),("Pericyte_SMC","ACTA2"),("Pericyte_SMC","TAGLN"),
]

marker_df = pd.DataFrame(markers, columns=["group","gene"])
marker_df["gene_upper"] = marker_df["gene"].str.upper()
marker_set = set(marker_df["gene_upper"])

print("[INFO] A_UMAP:", A_UMAP, flush=True)
print("[INFO] MATRIX:", MATRIX, flush=True)

umap = pd.read_csv(A_UMAP, sep="\t")
if "CB" not in umap.columns or "cluster" not in umap.columns:
    raise ValueError(f"UMAP columns={umap.columns.tolist()}, need CB and cluster")

umap["cluster"] = pd.to_numeric(umap["cluster"], errors="coerce").astype(int)
cell2cluster = dict(zip(umap["CB"].astype(str), umap["cluster"]))
cell_set = set(cell2cluster)

print("[INFO] selected cells:", len(cell_set), flush=True)
print("[INFO] cluster counts:", flush=True)
print(umap["cluster"].value_counts().sort_index().to_string(), flush=True)

# marker gene -> promoter region
bed = pd.read_csv(BED, sep="\t", header=None, dtype=str, comment="#")
records = []

for _, row in bed.iterrows():
    vals = [str(x) for x in row.tolist() if str(x) != "nan"]
    text = " ".join(vals)
    ens = None
    m = re.search(r"ENSSSCG[0-9]+", text)
    if m:
        ens = m.group(0)
    if ens is None:
        continue

    toks = set(t.upper() for t in re.split(r"[\s,;|:/()]+", text) if t)
    hits = marker_set.intersection(toks)
    for h in hits:
        for _, r in marker_df[marker_df["gene_upper"] == h].iterrows():
            records.append({
                "region_id": ens,
                "group": r["group"],
                "gene": r["gene"],
                "label": f'{r["group"]}:{r["gene"]}'
            })

region_map = pd.DataFrame(records).drop_duplicates()
region_map.to_csv(OUTDIR / "A_no01_marker_region_map.tsv", sep="\t", index=False)

print("[INFO] marker regions:", region_map.shape[0], flush=True)
print("[INFO] marker genes found:", region_map["gene"].nunique(), flush=True)

missing = sorted(marker_set - set(region_map["gene"].str.upper()))
pd.DataFrame({"missing_gene": missing}).to_csv(OUTDIR / "A_no01_missing_marker_genes.tsv", sep="\t", index=False)

if region_map.empty:
    raise RuntimeError("No marker region matched BED.")

region2label = region_map.groupby("region_id")["label"].apply(list).to_dict()
region_set = set(region2label)

# aggregate marker methylation by cluster
agg = {}
usecols = ["CB","region_id","modified_calls","covered_calls"]
reader = pd.read_csv(MATRIX, sep="\t", usecols=usecols, chunksize=1000000)

kept_total = 0
for idx, chunk in enumerate(reader, 1):
    chunk["CB"] = chunk["CB"].astype(str)
    chunk = chunk[chunk["CB"].isin(cell_set) & chunk["region_id"].isin(region_set)].copy()
    if chunk.empty:
        continue

    chunk["cluster"] = chunk["CB"].map(cell2cluster)
    chunk["modified_calls"] = pd.to_numeric(chunk["modified_calls"], errors="coerce").fillna(0)
    chunk["covered_calls"] = pd.to_numeric(chunk["covered_calls"], errors="coerce").fillna(0)

    rows = []
    for _, x in chunk.iterrows():
        for lab in region2label.get(x["region_id"], []):
            rows.append((int(x["cluster"]), lab, float(x["modified_calls"]), float(x["covered_calls"])))

    tmp = pd.DataFrame(rows, columns=["cluster","label","modified_calls","covered_calls"])
    g = tmp.groupby(["cluster","label"], as_index=False)[["modified_calls","covered_calls"]].sum()

    for _, r in g.iterrows():
        key = (int(r["cluster"]), r["label"])
        agg.setdefault(key, [0.0, 0.0])
        agg[key][0] += float(r["modified_calls"])
        agg[key][1] += float(r["covered_calls"])

    kept_total += len(chunk)
    print(f"[INFO] chunk {idx}, kept rows={len(chunk)}, total={kept_total}", flush=True)

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
    raise RuntimeError("No marker records found in matrix.")

# group-level
group_long = long.groupby(["cluster","group"], as_index=False)[["modified_calls","covered_calls"]].sum()
group_long["modified_rate"] = group_long["modified_calls"] / group_long["covered_calls"]
group_long["hypomethylation_activity"] = -group_long["modified_rate"]
group_long.to_csv(OUTDIR / "A_no01_marker_group_activity_long.tsv", sep="\t", index=False)

cluster_order = sorted(group_long["cluster"].unique())
group_order = [
    "Excitatory_neuron","Upper_layer_excitatory","Deep_layer_excitatory",
    "Inhibitory_neuron","Astrocyte","Oligodendrocyte","OPC",
    "Microglia","Endothelial","Pericyte_SMC"
]

gmat = group_long.pivot_table(index="group", columns="cluster", values="hypomethylation_activity", aggfunc="mean")
gmat = gmat.reindex(index=[g for g in group_order if g in gmat.index], columns=cluster_order)

arr = gmat.to_numpy(dtype=float)
mu = np.nanmean(arr, axis=1, keepdims=True)
sd = np.nanstd(arr, axis=1, keepdims=True)
sd[sd == 0] = np.nan
z = (arr - mu) / sd
z = np.nan_to_num(z, nan=0.0)
z = np.clip(z, -2, 2)

gzdf = pd.DataFrame(z, index=gmat.index, columns=[f"C{c}" for c in cluster_order])
gzdf.to_csv(OUTDIR / "A_no01_marker_group_activity_zscore.tsv", sep="\t")

# marker gene-level
mat = long.pivot_table(index="label", columns="cluster", values="hypomethylation_activity", aggfunc="mean")
label_order = []
for g, gene in markers:
    lab = f"{g}:{gene}"
    if lab in mat.index:
        label_order.append(lab)
mat = mat.reindex(index=label_order, columns=cluster_order)

arr = mat.to_numpy(dtype=float)
mu = np.nanmean(arr, axis=1, keepdims=True)
sd = np.nanstd(arr, axis=1, keepdims=True)
sd[sd == 0] = np.nan
z = (arr - mu) / sd
z = np.nan_to_num(z, nan=0.0)
z = np.clip(z, -2, 2)

zdf = pd.DataFrame(z, index=mat.index, columns=[f"C{c}" for c in cluster_order])
zdf.to_csv(OUTDIR / "A_no01_marker_promoter_activity_zscore.tsv", sep="\t")

def plot_heatmap(df, outprefix, title, figsize):
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(df.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(df.shape[1]))
    ax.set_xticklabels(df.columns)
    ax.set_yticks(range(df.shape[0]))
    ax.set_yticklabels(df.index, fontsize=7 if df.shape[0] > 20 else 10)
    ax.set_title(title)
    ax.set_xlabel("A rerun cluster")
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Relative promoter hypomethylation activity\nrow z-score")
    plt.tight_layout()
    fig.savefig(str(outprefix) + ".png", dpi=300)
    fig.savefig(str(outprefix) + ".pdf")
    plt.close(fig)

plot_heatmap(
    gzdf,
    OUTDIR / "A_no01_marker_group_activity_heatmap",
    "A rerun UMAP: marker group activity",
    (7.2, 4.8)
)

plot_heatmap(
    zdf,
    OUTDIR / "A_no01_marker_promoter_activity_heatmap",
    "A rerun UMAP: promoter hypomethylation marker activity",
    (8, max(8, 0.23 * zdf.shape[0]))
)

top_rows = []
for c in gzdf.columns:
    s = gzdf[c].sort_values(ascending=False)
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
