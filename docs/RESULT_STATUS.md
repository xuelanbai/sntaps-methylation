# 结果状态与可靠性边界

## 当前可复现的最终结果

最终主图输入是去除旧 P7 cluster0/1 后重建的 `P7_no01` promoter residual UMAP：

- 7,958 个细胞进入最终 6 群结果。
- C0：Inhibitory neurons（1,069）
- C1：Excitatory neurons I（985）
- C2：Vascular-associated cells（1,448）
- C3：Oligodendrocyte precursor cells / 胶质相关趋势（1,210）
- C4：Astrocytes（1,252）
- C5：Excitatory neurons II（1,994）

最终 UMAP 与 reads_used、covered_calls、covered_bins、modified_rate 的 Spearman 相关绝对值均较低（文档记录约 0.006–0.056），说明最终二维坐标不像早期版本那样直接由单一 QC 指标拉开。

## 不作为正式主结果的路线

- 原始 100-kb bin UMAP。
- gene-body / gene-body-no-promoter 的 100-kb 变体。
- 仅按覆盖度优先的 observed-only promoter UMAP。
- 旧 P7 的 8 群结果本身。

这些路线要么和技术指标高度相关，要么矩阵缺失率超过 99%，要么强制聚类后没有稳定、可解释的分群。

## 仍未完成的验证

1. **缺失主线脚本**：构建 `P7_no01` promoter 矩阵的脚本未包含在本次代码包中，需从旧服务器 shell history、登录节点目录或结果日志恢复。
2. **参数敏感性**：文档明确记录尚未系统比较 cluster 数、PC 数、随机种子和邻居参数。
3. **去除旧 cluster0/1 的合理性**：已有 QC/marker 证据支持其为低信息主体，但仍需把保留/去除两版的结论并列验证。
4. **细胞类型注释**：当前以 promoter 低甲基化 marker 趋势为主要证据；C2/C3/C4 有混合 marker 信号，不能写成高置信度精细类型。
5. **环境版本**：旧服务器没有随代码包保存完整的显式版本锁定；迁移前要从旧环境导出或重新建立最小环境并跑小样本验收。

## 正式报告前的最低验收

- 同一输入、同一种子重跑最终 UMAP，细胞集合和 cluster 数完全一致。
- 至少用 3 个随机种子、2 组 PC 数、2 组邻居数检查主结构稳定性。
- 每个 cluster 同时检查 QC、技术组来源、marker heatmap 和 cluster 大小。
- 对去除旧 cluster0/1 前后结果做并列表和敏感性说明。
- 用一个小规模 lane 完成从 split 到 promoter matrix 的端到端 smoke test。
