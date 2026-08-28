#!/bin/bash
set -euo pipefail

cd /share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation

PY=/share/org/YZWL/yzwl_wujx/.conda/envs/snTAPS_analysis/bin/python

SCRIPT100=03_scripts/34_umap_100kb_qcregress_fast.py
SCRIPTREG=03_scripts/region_matrix_paper/58_umap_region_qcregress_fast.py

BIN=08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.bin_counts.tsv.gz
PROM=08_qc_summary/final_single_sample/region_matrix/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.region_counts.tsv.gz

QC2=08_qc_summary/final_single_sample/cell_filter_balanced_core/snTAPS-FL-1.balanced_core2000.global_per_cb.tsv
QC5=08_qc_summary/final_single_sample/cell_filter_balanced_core/snTAPS-FL-1.balanced_core5000.global_per_cb.tsv

OUTBASE=08_qc_summary/final_single_sample/region_umap_paper/filteredCB_promoter_100kb_rescue
LOGDIR=09_logs/region_umap_paper/filteredCB_promoter_100kb_rescue

mkdir -p "$OUTBASE" "$LOGDIR"

for f in "$PY" "$SCRIPT100" "$SCRIPTREG" "$BIN" "$PROM" "$QC2" "$QC5"; do
  if [ ! -e "$f" ]; then
    echo "[ERROR] missing: $f"
    exit 1
  fi
done

echo "[INFO] check qcregress args"
$PY "$SCRIPT100" -h | grep -q -- "--qc-keys" || {
  echo "[ERROR] $SCRIPT100 does not support --qc-keys"
  exit 2
}
$PY "$SCRIPTREG" -h | grep -q -- "--qc-keys" || {
  echo "[ERROR] $SCRIPTREG does not support --qc-keys"
  exit 2
}

submit_100kb () {
  local tag=$1
  local qc=$2
  local min_bin_cells=$3
  local top_bins=$4

  local name="all100kb_${tag}_seurat_qcreg"
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

$PY $SCRIPT100 \\
  --global-qc "$qc" \\
  --bin-counts "$BIN" \\
  --outdir "$outdir" \\
  --prefix "snTAPS-FL-1.${name}" \\
  --mode seurat_norm \\
  --qc-keys reads_used,modified_rate \\
  --min-reads 0 \\
  --min-covered-calls 0 \\
  --min-covered-bins 0 \\
  --cellbin-min-cov 1 \\
  --max-cells 200000 \\
  --min-bin-cells "$min_bin_cells" \\
  --top-variable-bins "$top_bins" \\
  --n-pcs 30 \\
  --n-clusters 8 \\
  --n-jobs 16 \\
  --seed 1
EOF

  echo "[INFO] submit 100kb $name"
  csub < "$shfile"
}

submit_promoter () {
  local tag=$1
  local qc=$2

  local name="promoter_${tag}_seurat_qcreg"
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

$PY $SCRIPTREG \\
  --region-counts "$PROM" \\
  --cell-qc "$qc" \\
  --outdir "$outdir" \\
  --prefix "snTAPS-FL-1.${name}" \\
  --mode seurat_norm \\
  --qc-keys reads_used,modified_rate \\
  --n-top-regions 1000 \\
  --min-region-cells 20 \\
  --min-regions-per-cell 3 \\
  --n-pcs 30 \\
  --n-clusters 8 \\
  --seed 1
EOF

  echo "[INFO] submit promoter $name"
  csub < "$shfile"
}

# 100kb：严格CB后再试一次
submit_100kb "core2000" "$QC2" 500 1000
submit_100kb "core5000" "$QC5" 1250 1000

# promoter：只做验证，不再大规模调参
submit_promoter "core2000" "$QC2"
submit_promoter "core5000" "$QC5"

echo "[INFO] submitted filteredCB promoter/100kb rescue jobs"
