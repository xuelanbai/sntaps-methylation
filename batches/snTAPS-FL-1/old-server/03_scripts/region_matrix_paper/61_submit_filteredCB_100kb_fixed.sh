#!/bin/bash
set -euo pipefail

cd /share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation

PY=/share/org/YZWL/yzwl_wujx/.conda/envs/snTAPS_analysis/bin/python
SCRIPT=03_scripts/34_umap_100kb_qcregress_fast.py

BIN=08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.bin_counts.tsv.gz

QC2=08_qc_summary/final_single_sample/cell_filter_balanced_core/snTAPS-FL-1.balanced_core2000.global_per_cb.tsv
QC5=08_qc_summary/final_single_sample/cell_filter_balanced_core/snTAPS-FL-1.balanced_core5000.global_per_cb.tsv

OUTBASE=08_qc_summary/final_single_sample/region_umap_paper/filteredCB_100kb_fixed
LOGDIR=09_logs/region_umap_paper/filteredCB_100kb_fixed
mkdir -p "$OUTBASE" "$LOGDIR"

for f in "$PY" "$SCRIPT" "$BIN" "$QC2" "$QC5"; do
  if [ ! -e "$f" ]; then
    echo "[ERROR] missing: $f"
    exit 1
  fi
done

submit_one () {
  local tag=$1
  local qc=$2
  local min_bin_cells=$3
  local top_bins=$4

  local name="all100kb_${tag}_seurat_qcreg_minbin${min_bin_cells}_top${top_bins}"
  local outdir="$OUTBASE/$name"
  local shfile="03_scripts/region_matrix_paper/run_${name}.csub"

  mkdir -p "$outdir"

  cat > "$shfile" <<EOF
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
  --global-qc "$qc" \\
  --bin-counts "$BIN" \\
  --outdir "$outdir" \\
  --prefix "snTAPS-FL-1.${name}" \\
  --mode seurat_norm \\
  --qc-keys reads_used,modified_rate \\
  --cellbin-min-cov 1 \\
  --min-bin-cells "$min_bin_cells" \\
  --top-variable-bins "$top_bins" \\
  --min-selected-bins 100 \\
  --n-pcs 30 \\
  --n-clusters 8 \\
  --n-neighbors 30 \\
  --min-dist 0.3 \\
  --n-components 2 \\
  --n-jobs 16 \\
  --seed 1
EOF

  echo "[INFO] submit $name"
  csub < "$shfile"
}

submit_one core2000 "$QC2" 500 1000
submit_one core5000 "$QC5" 1250 1000

echo "[INFO] submitted fixed 100kb jobs"
