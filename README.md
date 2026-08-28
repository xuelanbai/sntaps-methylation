# snTAPS methylation

猪额叶 snTAPS 单细胞/单核 TAPS 甲基化分析代码归档与迁移版。

本仓库同时保留两种内容：

- `workflow/`：以后迁移、换路径或重跑时优先使用的简化入口。
- `batches/snTAPS-FL-1/old-server/`：旧服务器原始脚本快照，用于追溯当时实际运行过的命令；其中包含成功主线、探索版本、失败测试和自动生成的批量作业，不能全部视为正式流程。
- `SOURCE_SNAPSHOT.tsv`：341 个归档文件的大小和 SHA-256，用于迁移后完整性核对。

## 当前推荐主线

```text
raw FASTQ
  -> barcode 拆分（T7 + T5 + I7 + I5；写入 CB/CN）
  -> fastp
  -> 按 lane 比对
  -> 按 CB 做 CpG/细胞级统计
  -> 每个技术组合并 8 个 lane
  -> 合并 8 个技术组
  -> strict2 细胞白名单
  -> promoter/TSS ±1.5 kb 矩阵
  -> 旧 P7 仅用于识别低信息 cluster0/1
  -> 去除旧 cluster0/1 后重建 P7_no01 promoter 矩阵
  -> residual + QC regression + dropPC1 + smooth100 UMAP
  -> promoter 低甲基化 marker 趋势注释
  -> 最终作图
```

详细脚本对应关系见 [docs/WORKFLOW.md](docs/WORKFLOW.md)，结果可用性和未完成项见 [docs/RESULT_STATUS.md](docs/RESULT_STATUS.md)。

## 重要结论

- 100-kb bin、gene-body 100-kb 等路线是探索过程，现有结果主要被覆盖度、测序深度或整体 modified rate 驱动，不作为最终主结果。
- 旧 P7 的 cluster0/1 占比异常大，后续主图使用去除它们后重建的 `P7_no01` 结果。
- 最终 6 群注释目前仍是“趋势性注释”，不是已经完全验证的精细细胞类型结论。
- 文档记录中用于构建 `P7_no01` promoter 矩阵的脚本没有出现在本次旧服务器代码包内；这是迁移前必须补齐的唯一主线代码缺口。

## 快速开始

```bash
cp workflow/config.example.env config.env
# 修改 config.env
bash workflow/run_final_analysis.sh preflight config.env
```

`run_final_analysis.sh` 只覆盖已经稳定下来的 promoter 矩阵、最终 UMAP、marker 和主图步骤。早期拆分、比对、calling 和合并的实际批次脚本仍在旧服务器快照中，待新服务器环境确定后再做第二轮参数化。

## 数据边界

仓库不包含 FASTQ、BAM、参考基因组、结果矩阵或个人环境目录。所有大文件都应保存在服务器项目目录，并通过 `config.env` 指向。
