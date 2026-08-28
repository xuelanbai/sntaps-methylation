import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--marker", required=True)
parser.add_argument("--agg", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--topn", type=int, default=6)
args = parser.parse_args()

marker_file = Path(args.marker)
agg_file = Path(args.agg)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

mk = pd.read_csv(marker_file, sep="\t")

# 选每个cluster最强的低甲基化promoter marker
mk = mk.copy()
mk["cluster"] = mk["cluster"].astype(str)

if "abs_delta" not in mk.columns:
    mk["abs_delta"] = mk["delta_rate"].abs()

def pick_label(row):
    for c in ["mapped_symbol", "gene_symbol"]:
        if c in row and pd.notna(row[c]):
            v = str(row[c])
            if v and not v.startswith("ENSSSCG"):
                return v
    rid = str(row["region_id"])
    if "|" in rid:
        return rid.split("|")[-1]
    return rid

mk["label"] = mk.apply(pick_label, axis=1)
mk = mk.sort_values(["cluster", "abs_delta"], ascending=[True, False])
sel = mk.groupby("cluster", group_keys=False).head(args.topn)

# 避免同一个region重复出现
sel = sel.drop_duplicates("region_id").copy()
selected_regions = set(sel["region_id"].astype(str))

sel.to_csv(outdir / "P7.selected_marker_regions_for_heatmap.tsv", sep="\t", index=False)

# 从cluster-region聚合表中只读取marker对应区域
chunks = []
for chunk in pd.read_csv(agg_file, sep="\t", chunksize=1000000):
    chunk["region_id"] = chunk["region_id"].astype(str)
    sub = chunk[chunk["region_id"].isin(selected_regions)].copy()
    if len(sub) > 0:
        chunks.append(sub)

if len(chunks) == 0:
    raise SystemExit("No selected marker regions found in agg table.")

agg = pd.concat(chunks, ignore_index=True)
agg["cluster"] = agg["cluster"].astype(str)
agg["modified_rate"] = agg["modified_calls"] / agg["covered_calls"].replace(0, np.nan)

rate = agg.pivot_table(
    index="region_id",
    columns="cluster",
    values="modified_rate",
    aggfunc="mean"
)

clusters = sorted(rate.columns, key=lambda x: int(x) if str(x).isdigit() else str(x))
rate = rate.reindex(columns=clusters)

# 按marker选择顺序排序
ordered_regions = [r for r in sel["region_id"].astype(str).tolist() if r in rate.index]
rate = rate.reindex(ordered_regions)

labels = sel.set_index("region_id").loc[ordered_regions, "label"].astype(str).tolist()
source_cluster = sel.set_index("region_id").loc[ordered_regions, "cluster"].astype(str).tolist()

# promoter低甲基化代表潜在活性高，所以用负的modified_rate做activity
activity = -rate

# 每个marker按cluster方向做z-score
mat = activity.copy()
mat = mat.sub(mat.mean(axis=1), axis=0)
mat = mat.div(mat.std(axis=1).replace(0, np.nan), axis=0)
mat = mat.fillna(0)

# 输出矩阵
rate.to_csv(outdir / "P7.marker_promoter_modified_rate_matrix.tsv", sep="\t")
mat.to_csv(outdir / "P7.marker_promoter_activity_zscore_matrix.tsv", sep="\t")

# 画热图
fig_h = max(5, 0.22 * len(labels))
fig, ax = plt.subplots(figsize=(6.5, fig_h))

im = ax.imshow(mat.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)

ax.set_xticks(np.arange(len(clusters)))
ax.set_xticklabels([f"Cluster {c}" for c in clusters], rotation=45, ha="right")

yticklabels = [f"{lab}  C{cl}" for lab, cl in zip(labels, source_cluster)]
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels(yticklabels, fontsize=7)

ax.set_title("Cluster-specific promoter hypomethylation markers")
ax.set_xlabel("Cell cluster")
ax.set_ylabel("Marker promoter")

cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
cbar.set_label("Relative promoter activity score\nz-score of -modified_rate")

plt.tight_layout()
plt.savefig(outdir / "P7.marker_promoter_activity_heatmap.png", dpi=300)
plt.savefig(outdir / "P7.marker_promoter_activity_heatmap.pdf")
plt.close()

print("[DONE]")
print(outdir / "P7.marker_promoter_activity_heatmap.png")
print(outdir / "P7.marker_promoter_activity_heatmap.pdf")
print(outdir / "P7.marker_promoter_activity_zscore_matrix.tsv")
