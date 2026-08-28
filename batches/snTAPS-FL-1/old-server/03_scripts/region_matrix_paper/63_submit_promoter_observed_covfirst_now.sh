#!/bin/bash
set -euo pipefail

cd /share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation

PY=/share/org/YZWL/yzwl_wujx/.conda/envs/snTAPS_analysis/bin/python
SCRIPT=03_scripts/region_matrix_paper/63_umap_promoter_observed_covfirst.py
REGION=08_qc_summary/final_single_sample/region_matrix/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.region_counts.tsv.gz

# 优先用最新expected_realcell；找不到再回退到balanced_core
QC2=$(find 08_qc_summary/final_single_sample -type f -name "*expected_realcell_top2000*global_per_cb.tsv" | head -1)
QC5=$(find 08_qc_summary/final_single_sample -type f -name "*expected_realcell_top5000*global_per_cb.tsv" | head -1)

[ -n "$QC2" ] || QC2=08_qc_summary/final_single_sample/cell_filter_balanced_core/snTAPS-FL-1.balanced_core2000.global_per_cb.tsv
[ -n "$QC5" ] || QC5=08_qc_summary/final_single_sample/cell_filter_balanced_core/snTAPS-FL-1.balanced_core5000.global_per_cb.tsv

OUTBASE=08_qc_summary/final_single_sample/region_umap_paper/promoter_observed_covfirst_now
LOGDIR=09_logs/region_umap_paper/promoter_observed_covfirst_now
mkdir -p "$OUTBASE" "$LOGDIR"

echo "[INFO] QC2=$QC2"
echo "[INFO] QC5=$QC5"

submit_one () {
  local tag=$1
  local qc=$2
  local min_region_cells=$3

  local name="promoter_${tag}_obs_covfirst_minCell${min_region_cells}_top1000_k4"
  local outdir="$OUTBASE/$name"
  local shfile="03_scripts/region_matrix_paper/run_${name}.csub"

  mkdir -p "$outdir"

  cat > "$shfile" <<EOF
#!/bin/bash
#CSUB -J ${name}
#CSUB -q c02
#CSUB -n 8
#CSUB -R 'span[hosts=1]'
#CSUB -cwd /share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation
#CSUB -o ${LOGDIR}/${name}.out
#CSUB -e ${LOGDIR}/${name}.err

set -euo pipefail
cd /share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation

$PY $SCRIPT \\
  --region-counts "$REGION" \\
  --cell-qc "$qc" \\
  --outdir "$outdir" \\
  --prefix "snTAPS-FL-1.${name}" \\
  --min-cov 1 \\
  --n-top-regions 1000 \\
  --min-region-cells "$min_region_cells" \\
  --min-regions-per-cell 3 \\
  --n-clusters 4 \\
  --seed 1
EOF

  echo "[INFO] submit $name"
  csub < "$shfile"
}

# 只跑三个，不再铺太多
submit_one core2000 "$QC2" 5
submit_one core5000 "$QC5" 5
submit_one core5000 "$QC5" 3

echo "[INFO] submitted promoter observed coverage-first jobs"
