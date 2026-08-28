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
MATRIX = Path("08_qc_summary/final_single_sample/region_matrix/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.region_counts.tsv.gz")
OUTDIR = Path("08_qc_summary/final_single_sample/marker_pig_frontal_lobe/A_no01_rerun_marker_v2_matrixID")
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
marker_df = pd.DataFrame(markers, columns=["group", "gene"])
marker_df["gene_upper"] = marker_df["gene"].str.upper()
marker_set = set(marker_df["gene_upper"])

print("[INFO] A_UMAP:", A_UMAP, flush=True)
print("[INFO] MATRIX:", MATRIX, flush=True)

# ---------- read A cells ----------
umap = pd.read_csv(A_UMAP, sep="\t")
if "CB" not in umap.columns or "cluster" not in umap.columns:
    raise ValueError(f"UMAP columns={umap.columns.tolist()}, need CB and cluster")

umap["CB"] = umap["CB"].astype(str)
umap["cluster"] = pd.to_numeric(umap["cluster"], errors="coerce").astype(int)
cell2cluster = dict(zip(umap["CB"], umap["cluster"]))
cell_set = set(cell2cluster)

print("[INFO] selected cells:", len(cell_set), flush=True)
print("[INFO] cluster counts:", flush=True)
print(umap["cluster"].value_counts().sort_index().to_string(), flush=True)

# ---------- collect matrix region IDs ----------
print("[INFO] collect matrix region IDs", flush=True)
matrix_regions = set()
for i, chunk in enumerate(pd.read_csv(MATRIX, sep="\t", usecols=["region_id"], chunksize=1000000), 1):
    matrix_regions.update(chunk["region_id"].astype(str).unique())
    if i % 10 == 0:
        print(f"[INFO] scanned chunks={i}, unique matrix regions={len(matrix_regions)}", flush=True)

print("[INFO] total matrix regions:", len(matrix_regions), flush=True)

# ---------- parse BED, but choose only IDs existing in matrix ----------
print("[INFO] parse BED with matrix-compatible region_id", flush=True)
records = []
with open(BED) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = re.split(r"\s+", line.strip())
        text = " ".join(parts)

        ens_all = re.findall(r"ENSSSCG[0-9]+", text)
        ens_in_matrix = [x for x in ens_all if x in matrix_regions]
        if not ens_in_matrix:
            continue

        toks = set(t.upper() for t in re.split(r"[\s,;|:/()]+", text) if t)
        hit_genes = marker_set.intersection(toks)
        if not hit_genes:
            continue

        # 一般每行只应该有一个真正region_id在矩阵中；如果多个，全部保留
        for rid in sorted(set(ens_in_matrix)):
            for g in sorted(hit_genes):
                sub = marker_df[marker_df["gene_upper"] == g]
                for _, r in sub.iterrows():
                    records.append({
                        "region_id": rid,
                        "group": r["group"],
                        "gene": r["gene"],
                        "label": f'{r["group"]}:{r["gene"]}'
                    })

region_map = pd.DataFrame(records).drop_duplicates()
region_map.to_csv(OUTDIR / "A_no01_marker_region_map.matrixID.tsv", sep="\t", index=False)

found = set(region_map["gene"].str.upper()) if not region_map.empty else set()
missing = sorted(marker_set - found)
pd.DataFrame({"missing_gene": missing}).to_csv(OUTDIR / "A_no01_missing_marker_genes.matrixID.tsv", sep="\t", index=False)

print("[INFO] marker regions compatible with matrix:", region_map["region_id"].nunique() if not region_map.empty else 0, flush=True)
print("[INFO] marker genes found:", region_map["gene"].nunique() if not region_map.empty else 0, flush=True)
print("[INFO] missing marker genes:", len(missing), flush=True)

if region_map.empty:
    raise RuntimeError("No marker regions matched matrix region IDs. Need inspect BED columns.")

region2label = region_map.groupby("region_id")["label"].apply(list).to_dict()
region_set = set(region2label)

# ---------- stream matrix and aggregate ----------
print("[INFO] scan matrix for A cells x marker regions", flush=True)
agg = {}
usecols = ["CB", "region_id", "modified_calls", "covered_calls"]

total_cb_hit = 0
total_region_hit = 0
total_both = 0

for i, chunk in enumerate(pd.read_csv(MATRIX, sep="\t", usecols=usecols, chunksize=1000000), 1):
    chunk["CB"] = chunk["CB"].astype(str)
    chunk["region_id"] = chunk["region_id"].astype(str)

    cb_mask = chunk["CB"].isin(cell_set)
    rg_mask = chunk["region_id"].isin(region_set)
    total_cb_hit += int(cb_mask.sum())
    total_region_hit += int(rg_mask.sum())

    sub = chunk[cb_mask & rg_mask].copy()
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
        for lab in region2label.get(x["region_id"], []):
            rows.append((int(x["cluster"]), lab, float(x["modified_calls"]), float(x["covered_calls"])))

    tmp = pd.DataFrame(rows, columns=["cluster", "label", "modified_calls", "covered_calls"])
    g = tmp.groupby(["cluster", "label"], as_index=False)[["modified_calls", "covered_calls"]].sum()

    for _, r in g.iterrows():
        key = (int(r["cluster"]), r["label"])
        agg.setdefault(key, [0.0, 0.0])
        agg[key][0] += float(r["modified_calls"])
        agg[key][1] += float(r["covered_calls"])

    print(f"[INFO] chunk={i}, kept={len(sub)}, total_both={total_both}", flush=True)

print("[INFO] total rows with A cells:", total_cb_hit, flush=True)
print("[INFO] total rows with marker regions:", total_region_hit, flush=True)
print("[INFO] total rows with A cells AND marker regions:", total_both, flush=True)

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
    raise RuntimeError("No marker records found after matrix-compatible ID mapping.")

# ---------- group-level ----------
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
    z = np.nan_to_num(z, nan=0.0)
    z = np.clip(z, -2, 2)
    return pd.DataFrame(z, index=mat.index, columns=[f"C{c}" for c in mat.columns])

gmat = group_long.pivot_table(index="group", columns="cluster", values="hypomethylation_activity", aggfunc="mean")
gmat = gmat.reindex(index=[g for g in group_order if g in gmat.index], columns=cluster_order)
gzdf = row_zscore(gmat)
gzdf.to_csv(OUTDIR / "A_no01_marker_group_activity_zscore.tsv", sep="\t")

# ---------- gene-level ----------
mat = long.pivot_table(index="label", columns="cluster", values="hypomethylation_activity", aggfunc="mean")
label_order = []
for g, gene in markers:
    lab = f"{g}:{gene}"
    if lab in mat.index:
        label_order.append(lab)
mat = mat.reindex(index=label_order, columns=cluster_order)
zdf = row_zscore(mat)
zdf.to_csv(OUTDIR / "A_no01_marker_promoter_activity_zscore.tsv", sep="\t")

# ---------- delta summary ----------
delta_rows = []
for group, row in gmat.iterrows():
    s = row.dropna()
    if len(s) == 0:
        continue
    delta_rows.append({
        "group": group,
        "best_cluster": int(s.idxmax()),
        "max_activity": s.max(),
        "min_activity": s.min(),
        "delta": s.max() - s.min()
    })
pd.DataFrame(delta_rows).sort_values("delta", ascending=False).to_csv(
    OUTDIR / "A_no01_marker_group_delta_summary.tsv", sep="\t", index=False
)

# ---------- plots ----------
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

# ---------- top marker groups ----------
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
