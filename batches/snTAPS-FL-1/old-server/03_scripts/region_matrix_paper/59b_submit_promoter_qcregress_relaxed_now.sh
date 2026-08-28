#!/bin/bash
set -euo pipefail

cd /share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation

PY=/share/org/YZWL/yzwl_wujx/.conda/envs/snTAPS_analysis/bin/python
SCRIPT=03_scripts/region_matrix_paper/58_umap_region_qcregress_fast.py

# 这次不用cov10，用原始promoter region_counts，避免区域被筛到只剩20个
PROM=08_qc_summary/final_single_sample/region_matrix/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.region_counts.tsv.gz

QC=08_qc_summary/final_single_sample/cell_filter_strict_qc_core/snTAPS-FL-1.strict_qc_core_q25_q75_ge100k.core2000.global_per_cb.tsv

OUTBASE=08_qc_summary/final_single_sample/region_umap_paper/promoter_qcregress_rescue_relaxed
LOGDIR=09_logs/region_umap_paper/promoter_qcregress_rescue_relaxed
mkdir -p "$OUTBASE" "$LOGDIR"

submit_one () {
  local name=$1
  local mode=$2
  local qckeys=$3

  local outdir="$OUTBASE/$name"
  mkdir -p "$outdir"

  local sh="03_scripts/region_matrix_paper/run_${name}.csub"

  cat > "$sh" <<EOF
#!/bin/bash
#CSUB -J ${name}
#CSUB -q c02
#CSUB -n 16
#CSUB -R 'span[hosts=1]'
#CSUB -cwd /share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation
#CSUB -o ${LOGDIR}/${name}.out
#CSUB -e ${LOGDIR}/${name}.err

set -euo pipefail
cd /share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation

$PY $SCRIPT \\
  --region-counts "$PROM" \\
  --cell-qc "$QC" \\
  --outdir "$outdir" \\
  --prefix "snTAPS-FL-1.${name}" \\
  --mode "$mode" \\
  --qc-keys "$qckeys" \\
  --n-top-regions 3000 \\
  --min-region-cells 20 \\
  --min-regions-per-cell 5 \\
  --n-pcs 30 \\
  --n-clusters 8 \\
  --seed 1
EOF

  csub < "$sh"
}

submit_one "promoter_raw_q25q75_core2000_residual_depthOnly_relaxed" "residual" "reads_used"
submit_one "promoter_raw_q25q75_core2000_seurat_depthMod_relaxed" "seurat_norm" "reads_used,modified_rate"
