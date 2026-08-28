# 流程与脚本地图

## A. 已确认的主线

| 阶段 | 旧服务器脚本 | 主要输入 | 主要输出 | 状态 |
|---|---|---|---|---|
| 原始 FASTQ 链接/清单 | `03_scripts/00_link_raw_fastq.sh`, `01_make_fastq_lists.sh` | 测序公司释放目录 | `01_rawdata/`, FASTQ 清单 | 主线 |
| barcode 结构确认 | `07_extract_read_head_check.py`, `08_check_sntaps_barcode_structure.py` | 原始 R1/R2/header index | T5/T7/I5/I7 匹配统计 | 主线前检查 |
| barcode 拆分 | `09_split_sntaps_fastq.py`, `10_make_split_pair_table.sh`, `11_split_sntaps_all.csub` | 64 对 lane FASTQ、4 个 index 表 | split FASTQ、CB/CN header、拆分率 | 主线 |
| split 后质控 | `12_make_fastp_after_split_pair_table.sh`, `13_fastp_after_split.64core.csub`, `14_fastqc_after_split_fastp.csub` | split FASTQ | fastp FASTQ/HTML/JSON、FastQC | 主线 |
| lane 比对 | `16_make_align_pair_table.sh`, `17_bwa_align_by_lane.64core.cleanheader.csub`, `18_clean_fastq_header_for_bwaC.py` | fastp FASTQ、猪参考基因组 | 按 lane 排序 BAM/BAI | 主线 |
| TAPS calling | `20_call_taps_cpg_bins_by_cb.py`, `call_taps_full_by_lane/*.csub` | BAM、FASTA、CpG bins、CB 白名单 | lane 级 per-CB/bins/summary | 主线；64 个批次作业 |
| lane 合并 | `25_merge_lane_taps_one_sample.sh`, `merge_lane_results/*.csub` | 每技术组 8 lane | 技术组级 per-CB/bins/summary | 主线 |
| 技术组合并 | `final_merge/27_merge_all_technical_groups_one_sample.csub` | 8 个技术组 | `all_technical_merged.*` | 主线 |
| strict2 细胞筛选 | `make_strict_qc_core_cells.py` 及相关 QC 脚本 | 合并后的 cell-level QC | `strict2_CB.whitelist.txt`, `strict2_cells.global_qc.tsv` | 主线；阈值见下 |
| promoter 区域 | `reference/02_make_pig_gene_regions.py` | 猪 GTF/染色体长度 | `pig_promoter_TSS_1.5kb.bed` | 主线 |
| promoter 矩阵 | `region_matrix/40_call_region_matrix_by_cb.py`, `40_call_promoter_matrix_strict2.csub` | 64 BAM、promoter BED、strict2 白名单 | cell × promoter counts | 主线 |
| 旧 P7 | `region_matrix_paper/44_umap_region_matrix_paper_fast.py`, `46_P7_*.csub` | strict2 promoter 矩阵、global QC | 旧 P7 UMAP/cluster | 仅用于识别待排除群 |
| P7_no01 矩阵 | 文档记录为 `region_matrix_paper/99_make_P7_no01_promoter_matrix.py` | 旧 P7 cluster、strict2 promoter 矩阵 | 去除旧 cluster0/1 后的矩阵 | **代码包缺失，必须补齐** |
| 最终 UMAP | `region_matrix_paper/100_P7_no01_rerun_promoter_residual_dropPC1_smooth100.csub` + `44_umap_region_matrix_paper_fast.py` | P7_no01 矩阵、global QC | 最终 6 群 UMAP | 主线 |
| marker 趋势 | `annotation/101_A_no01_rerun_marker_heatmap.py` | 最终 UMAP、promoter 矩阵、marker 表 | marker activity / Z-score / heatmap | 主线；注释仍需验证 |
| 主图 | `figures/110_make_A_main_figures.py` | 最终 UMAP、marker Z-score | Fig1–Fig3 | 主线 |

所有路径均相对于：

```text
batches/snTAPS-FL-1/old-server/
```

## B. 关键数据结构

- F1–F8 是同一样本的 8 个技术组；每个技术组有 L1–L8。
- read 前端结构：R1 前 8 bp = T5，R2 前 8 bp = T7，header index1 = I7，index2 = I5。
- 完整 barcode：`T7 + T5 + I7 + I5`，写入 FASTQ header 的 `CB:Z:`；拆分后去掉 R1/R2 左端 30 bp。
- 每个技术组先合并 8 个 lane，再将 8 个技术组合为最终样本 `snTAPS-FL-1`。

## C. 关键阈值与最终 UMAP 参数

promoter 矩阵：

- MAPQ ≥ 30
- base quality ≥ 20
- 每个 cell-promoter 至少 5 个 covered calls

最终 P7_no01 UMAP：

- `mode=residual`
- `max_cells=0`
- `min_cell_regions=500`
- `min_region_cells=1000`
- `top_variable_regions=2000`
- QC regression
- drop PC1
- 10 PCs
- 100 neighbors
- `min_dist=0.4`
- 6 clusters

## D. 探索/历史目录

- `region_matrix_paper/` 中 P1–P13、100-kb/gene-body、observed-only、balanced-core、qcregress A/B 等是参数和特征探索，不应和最终 A 版主线混用。
- `archive/wrong_order_fastp_before_split/` 是先 fastp 后拆分的错误顺序，正式分析禁止使用。
- `archive/failed_test_scripts/` 是已被 Python 脚本替代的早期测试。
- `archive/reference_scATAC_scripts/` 只用于来源追溯。
- `call_taps_full_by_lane/` 和 `merge_lane_results/` 是本批次自动展开的作业，保留用于复现，不是通用模板。
