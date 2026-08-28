# 旧服务器环境观察

代码中出现四个环境：

| 环境 | 观察到的用途 | 主要工具/包 |
|---|---|---|
| `snTAPS_qc` | 原始/拆分后 QC | fastp, FastQC, MultiQC |
| `snTAPS_align` | 比对和 BAM 处理 | bwa, samtools, sambamba/picard（部分脚本） |
| `snTAPS_analysis` | TAPS calling、矩阵、UMAP、作图 | pysam, pyfaidx, pandas, numpy, scipy, scikit-learn, scanpy, anndata, umap-learn, matplotlib, seaborn, openpyxl, pyarrow, h5py |
| `foodie_env` | 个别历史参考、bedGraphToBigWig 等工具 | 不应继续作为 snTAPS 主环境依赖 |

另外观察到：bedtools、mosdepth、deeptools、seqkit、pigz/parallel、chromap 等工具名，但并非全部属于最终主线。

## 迁移原则

- 新服务器至少拆成 QC、alignment、analysis 三个环境，不再从 `foodie_env` 借工具。
- 先用旧服务器导出的显式版本建立 lock 文件；如果无法导出，再按主线脚本的最小依赖重建。
- 所有环境安装后记录 `conda list --explicit`、`conda env export --from-history` 和关键工具 `--version`。
- 不要直接把旧服务器 51 处绝对 Python 路径复制到新环境；统一通过激活环境后的 `python3` 调用。
