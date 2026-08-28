# 迁移检查清单

- [ ] 从旧服务器恢复 `99_make_P7_no01_promoter_matrix.py`（或等价命令）。
- [ ] 导出 `snTAPS_qc`、`snTAPS_align`、`snTAPS_analysis` 的版本信息。
- [ ] 确认新服务器 CSUB 队列、CPU、内存和并行方式。
- [ ] 准备猪参考 FASTA、GTF、BWA index、FAI、chrom sizes、CpG bins、promoter BED。
- [ ] 核对 T5/T7/I5/I7 字典及其来源版本。
- [ ] 用单 lane 小样本验证 barcode 拆分率、read 长度和 CB header。
- [ ] 用单 lane 验证 BAM 中 CB 标签、比对率和 TAPS calling。
- [ ] 用一个技术组验证 8 lane 合并。
- [ ] 重跑 strict2、promoter matrix 和最终 P7_no01。
- [ ] 做参数敏感性和注释复核。
