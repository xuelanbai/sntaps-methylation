#!/bin/bash
set -euo pipefail

cd /share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation

PY=/share/org/YZWL/yzwl_wujx/.conda/envs/snTAPS_analysis/bin/python
SCRIPT=03_scripts/region_matrix_paper/58_umap_region_qcregress_fast.py

PROM=08_qc_summary/final_single_sample/region_matrix_covfilter/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.cov10.region_counts.tsv.gz

QC_Q20=08_qc_summary/final_single_sample/cell_filter_strict_qc_core/snTAPS-FL-1.strict_qc_core_q20_q80_ge100k.core2000.global_per_cb.tsv
QC_Q25=08_qc_summary/final_single_sample/cell_filter_strict_qc_core/snTAPS-FL-1.strict_qc_core_q25_q75_ge100k.core2000.global_per_cb.tsv

OUTBASE=08_qc_summary/final_single_sample/region_umap_paper/promoter_qcregress_rescue
LOGDIR=09_logs/region_umap_paper/promoter_qcregress_rescue
mkdir -p "$OUTBASE" "$LOGDIR"

submit_one () {
  local name=$1
  local qc=$2
  local mode=$3
  local qckeys=$4

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
  --cell-qc "$qc" \\
  --outdir "$outdir" \\
  --prefix "snTAPS-FL-1.${name}" \\
  --mode "$mode" \\
  --qc-keys "$qckeys" \\
  --n-top-regions 3000 \\
  --min-region-cells 100 \\
  --min-regions-per-cell 20 \\
  --n-pcs 30 \\
  --n-clusters 8 \\
  --seed 1
EOF

  csub < "$sh"
}

submit_one "promoter_cov10_q25q75_core2000_residual_depthOnly" "$QC_Q25" "residual" "reads_used"
submit_one "promoter_cov10_q25q75_core2000_seurat_depthMod" "$QC_Q25" "seurat_norm" "reads_used,modified_rate"
submit_one "promoter_cov10_q20q80_core2000_residual_depthOnly" "$QC_Q20" "residual" "reads_used"
submit_one "promoter_cov10_q20q80_core2000_seurat_depthMod" "$QC_Q20" "seurat_norm" "reads_used,modified_rate"
