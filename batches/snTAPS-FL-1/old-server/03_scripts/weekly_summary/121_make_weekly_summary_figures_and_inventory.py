#!/usr/bin/env python3
import os
import glob
import subprocess
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


PROJECT = Path("/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation")
os.chdir(PROJECT)

OUT = Path("09_results/weekly_summary_snTAPS_pig_frontal")
OUT.mkdir(parents=True, exist_ok=True)

# ===== 固定使用的最终文件，不要随便改 =====
OLD_UMAP = Path(
    "08_qc_summary/final_single_sample/region_umap_paper/"
    "P7_promoter_residual_all_dropPC1_smooth100/"
    "snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"
)

A_UMAP = Path(
    "08_qc_summary/final_single_sample/region_umap_paper/"
    "P7_no01_promoter_residual_dropPC1_smooth100/"
    "snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"
)

GLOBAL_QC = Path(
    "08_qc_summary/final_single_sample/"
    "snTAPS-FL-1.all_technical_merged.global_per_cb.tsv"
)

REGION_MATRIX_COV10 = Path(
    "08_qc_summary/final_single_sample/region_matrix_covfilter/"
    "snTAPS-FL-1.strict2.promoter_TSS_1.5kb.cov10.region_counts.tsv.gz"
)

REGION_MATRIX_A = Path(
    "08_qc_summary/final_single_sample/region_matrix_P7_no01/"
    "snTAPS-FL-1.P7_no01.promoter_TSS_1.5kb.region_counts.tsv.gz"
)

SITE_CONVERSION = Path(
    "08_qc_summary/final_single_sample/site_conversion_merged/"
    "snTAPS-FL-1.all_chr.CpG_site_conversion.tsv.gz"
)

CLUSTER_BW_DIR = Path(
    "08_qc_summary/final_single_sample/cluster_bigwig_A_no01_from_bam_ref/merged"
)

MARKER_HEATMAP_DIR = Path(
    "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/A_no01_rerun_marker"
)

# ===== 细胞类型注释：A版最终结果 =====
CELLTYPE = {
    0: "Inhibitory neurons",
    1: "Excitatory neurons I",
    2: "Vascular-associated cells",
    3: "Oligodendrocyte precursor cells",
    4: "Astrocytes",
    5: "Excitatory neurons II",
}

# 为了图上排版，长名字换行
CELLTYPE_SHORT = {
    0: "Inhibitory\nneurons",
    1: "Excitatory\nneurons I",
    2: "Vascular-associated\ncells",
    3: "Oligodendrocyte\nprecursor cells",
    4: "Astrocytes",
    5: "Excitatory\nneurons II",
}


def human_size(path):
    path = Path(path)
    if not path.exists():
        return "NA"
    try:
        out = subprocess.check_output(["du", "-sh", str(path)], text=True).split()[0]
        return out
    except Exception:
        return str(path.stat().st_size)


def exists(path):
    if path is None or str(path) == "":
        return "NA"
    return "yes" if Path(path).exists() else "no"


def first_existing(patterns):
    for p in patterns:
        hits = glob.glob(p)
        if hits:
            return hits[0]
    return ""


def add_inventory(rows, step, purpose, input_file, output_file, code_file, parameter, note):
    rows.append({
        "step": step,
        "purpose": purpose,
        "input_file_or_dir": str(input_file),
        "input_exists": exists(input_file) if input_file else "NA",
        "input_size": human_size(input_file) if input_file and Path(input_file).exists() else "NA",
        "output_file_or_dir": str(output_file),
        "output_exists": exists(output_file) if output_file else "NA",
        "output_size": human_size(output_file) if output_file and Path(output_file).exists() else "NA",
        "code_file": str(code_file),
        "code_exists": exists(code_file) if code_file else "NA",
        "key_parameter": parameter,
        "note": note,
    })


# ===== 0. 流程文件清单表 =====
rows = []

RAW_DIR = first_existing([
    "01.RawData/snTAPS-FL-1_snTAPS-FL-1-*",
    "01_rawdata/snTAPS-FL-1_snTAPS-FL-1-*",
    "01_rawdata/*snTAPS*",
    "01.RawData/*snTAPS*",
])

BAM_DIR = Path("06_align/bam_by_lane")

add_inventory(
    rows,
    "01_raw_fastq",
    "原始FASTQ输入，测序公司释放后的正式样本数据",
    RAW_DIR,
    RAW_DIR,
    "",
    "raw paired FASTQ by lane",
    "只记录正式样本目录，不把Undetermined作为主分析输入",
)

add_inventory(
    rows,
    "02_alignment",
    "按lane比对得到sorted BAM",
    RAW_DIR,
    BAM_DIR,
    "03_scripts/call_taps_full_by_lane",
    "lane-wise mapping; BAM contains CB tag",
    "后续cluster-level bw从BAM重新统计CpG信号时使用该目录",
)

add_inventory(
    rows,
    "03_site_conversion",
    "全样本CpG位点甲基化/转化统计",
    BAM_DIR,
    SITE_CONVERSION,
    "03_scripts/bigwig_1bp",
    "CpG-level modified_calls/covered_calls/modified_rate",
    "整体CpG位点统计结果，后续bigWig和质控均以此类结果为基础",
)

add_inventory(
    rows,
    "04_global_cell_QC",
    "单细胞全局QC统计",
    SITE_CONVERSION,
    GLOBAL_QC,
    "03_scripts",
    "per-CB reads_used/covered_calls/covered_bins/modified_rate",
    "用于UMAP前细胞过滤和后续QC图",
)

add_inventory(
    rows,
    "05_promoter_region_matrix_cov10",
    "promoter/TSS±1.5kb区域矩阵",
    GLOBAL_QC,
    REGION_MATRIX_COV10,
    "03_scripts/region_matrix_paper",
    "promoter_TSS_1.5kb; cov10",
    "旧P7和cluster0/1质控使用的主要region matrix",
)

add_inventory(
    rows,
    "06_old_P7_UMAP",
    "旧P7 UMAP，用于识别cluster0/1低信息主体",
    REGION_MATRIX_COV10,
    OLD_UMAP,
    "03_scripts/region_matrix_paper/44_umap_region_matrix_paper_fast.py",
    "promoter residual; dropPC1; smooth100",
    "仅用于解释为什么去除旧cluster0/1，不作为最终主图",
)

add_inventory(
    rows,
    "07_A_rerun_region_matrix",
    "去除旧cluster0/1后重新构建的promoter矩阵",
    OLD_UMAP,
    REGION_MATRIX_A,
    "03_scripts/region_matrix_paper",
    "remove old cluster0/1; promoter_TSS_1.5kb",
    "A版最终UMAP上游输入之一",
)

add_inventory(
    rows,
    "08_A_rerun_UMAP",
    "最终主图使用的A版UMAP",
    REGION_MATRIX_A,
    A_UMAP,
    "03_scripts/region_matrix_paper/44_umap_region_matrix_paper_fast.py",
    "promoter residual; dropPC1; smooth100; n_clusters=6",
    "本周总结主图只使用该UMAP",
)

add_inventory(
    rows,
    "09_marker_heatmap",
    "猪额叶marker趋势注释",
    A_UMAP,
    MARKER_HEATMAP_DIR,
    "03_scripts/annotation/101_A_no01_rerun_marker_heatmap.py",
    "promoter hypomethylation marker activity",
    "marker gene验证尚未完成，暂不作为本周主图",
)

add_inventory(
    rows,
    "10_cluster_bigwig",
    "每个A版cluster合并生成methyl/hypo/coverage bw",
    A_UMAP,
    CLUSTER_BW_DIR,
    "03_scripts/bigwig_1bp/113_make_cluster_cpg_bg_from_bam_ref_onechr.py; 03_scripts/bigwig_1bp/114_merge_cluster_bedGraph_to_bw.csub",
    "CB tag + reference + read sequence; mincov5; canonical chromosomes",
    "IGV辅助验证使用，注意解释时必须结合coverage",
)

inventory = pd.DataFrame(rows)
inventory.to_csv(OUT / "00_workflow_final_input_output_inventory.tsv", sep="\t", index=False)


# ===== 1. 读入最终A版UMAP =====
if not A_UMAP.exists():
    raise FileNotFoundError(f"Missing final A UMAP: {A_UMAP}")

df = pd.read_csv(A_UMAP, sep="\t")
df["cluster"] = df["cluster"].astype(int)
df["cluster_label"] = df["cluster"].map(lambda x: f"C{x}")
df["cell_type"] = df["cluster"].map(CELLTYPE)
df["cell_type_short"] = df["cluster"].map(CELLTYPE_SHORT)

counts = (
    df.groupby(["cluster", "cluster_label", "cell_type"])
    .size()
    .reset_index(name="n_cells")
    .sort_values("cluster")
)
counts["fraction"] = counts["n_cells"] / counts["n_cells"].sum()
counts.to_csv(OUT / "final_cluster_cell_counts.tsv", sep="\t", index=False)


# ===== 图1：最终UMAP+cluster数量 =====
fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), gridspec_kw={"width_ratios": [2.2, 1]})

ax = axes[0]
for c in sorted(df["cluster"].unique()):
    sub = df[df["cluster"] == c]
    ax.scatter(sub["UMAP1"], sub["UMAP2"], s=3, linewidths=0, label=CELLTYPE[c], alpha=0.9)

# 用中位数标注
for c in sorted(df["cluster"].unique()):
    sub = df[df["cluster"] == c]
    ax.text(
        sub["UMAP1"].median(),
        sub["UMAP2"].median(),
        CELLTYPE_SHORT[c],
        ha="center",
        va="center",
        fontsize=9,
        weight="bold",
        color="black",
    )

ax.set_title("Porcine frontal lobe single-cell methylation UMAP", fontsize=14)
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax = axes[1]
ax.barh(counts["cell_type"], counts["n_cells"])
ax.invert_yaxis()
ax.set_xlabel("Number of cells")
ax.set_title("Cluster cell counts")
for i, v in enumerate(counts["n_cells"]):
    ax.text(v, i, f" {v}", va="center", fontsize=8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(OUT / "Fig1_final_umap_and_cluster_counts.png", dpi=300)
plt.savefig(OUT / "Fig1_final_umap_and_cluster_counts.pdf")
plt.close()


# ===== 图2：最终6个cluster的QC =====
qc_cols = [x for x in ["reads_used", "covered_calls", "covered_bins", "modified_rate"] if x in df.columns]

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes = axes.flatten()

for ax, col in zip(axes, qc_cols):
    data = [df.loc[df["cluster"] == c, col].dropna().values for c in sorted(df["cluster"].unique())]
    labels = [f"C{c}" for c in sorted(df["cluster"].unique())]
    ax.boxplot(data, labels=labels, showfliers=False)
    ax.set_title(col)
    ax.set_xlabel("A rerun cluster")
    ax.set_ylabel(col)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

for ax in axes[len(qc_cols):]:
    ax.axis("off")

plt.suptitle("QC metrics across final methylation clusters", fontsize=14)
plt.tight_layout()
plt.savefig(OUT / "Fig2_final_cluster_QC_panel.png", dpi=300)
plt.savefig(OUT / "Fig2_final_cluster_QC_panel.pdf")
plt.close()


# ===== 图3：旧cluster0/1质控 =====
if OLD_UMAP.exists():
    old = pd.read_csv(OLD_UMAP, sep="\t")
    old["cluster"] = old["cluster"].astype(int)
    old["group"] = old["cluster"].map(lambda x: "old C0/C1\nlow-information major population" if x in [0, 1] else "old C2-C7\nretained island population")

    old_counts = old.groupby("cluster").size().reset_index(name="n_cells").sort_values("cluster")
    old_counts.to_csv(OUT / "old_P7_cluster_counts.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    ax = axes[0, 0]
    # old UMAP：0/1灰色，2-7彩色
    for c in sorted(old["cluster"].unique()):
        sub = old[old["cluster"] == c]
        alpha = 0.35 if c in [0, 1] else 0.9
        label = f"C{c}"
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=2, linewidths=0, alpha=alpha, label=label)
    ax.set_title("Old P7 UMAP colored by cluster")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(markerscale=3, fontsize=7, frameon=False, ncol=2)

    ax = axes[0, 1]
    ax.bar(old_counts["cluster"].astype(str), old_counts["n_cells"])
    ax.set_title("Old P7 cluster cell counts")
    ax.set_xlabel("Old P7 cluster")
    ax.set_ylabel("Number of cells")
    for i, v in enumerate(old_counts["n_cells"]):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=7, rotation=90)

    ax = axes[0, 2]
    group_counts = old.groupby("group").size().reset_index(name="n_cells")
    ax.bar(group_counts["group"], group_counts["n_cells"])
    ax.set_title("Old C0/C1 vs retained islands")
    ax.set_ylabel("Number of cells")
    ax.tick_params(axis="x", labelrotation=15)

    old_qc_cols = [x for x in ["reads_used", "covered_calls", "covered_bins"] if x in old.columns]
    for ax, col in zip(axes[1, :], old_qc_cols):
        data = [
            old.loc[old["cluster"].isin([0, 1]), col].dropna().values,
            old.loc[~old["cluster"].isin([0, 1]), col].dropna().values,
        ]
        ax.boxplot(data, labels=["old C0/C1", "old C2-C7"], showfliers=False)
        ax.set_title(col)
        ax.set_ylabel(col)

    for ax in axes.flatten():
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.suptitle("QC evidence for excluding old cluster0/1", fontsize=14)
    plt.tight_layout()
    plt.savefig(OUT / "Fig3_old_cluster01_QC_panel.png", dpi=300)
    plt.savefig(OUT / "Fig3_old_cluster01_QC_panel.pdf")
    plt.close()
else:
    print(f"[WARN] OLD_UMAP not found: {OLD_UMAP}")


print("[DONE] wrote:", OUT)
print("\n[cluster counts]")
print(counts)
print("\n[inventory]")
print(inventory[["step", "output_file_or_dir", "output_exists", "output_size"]])
