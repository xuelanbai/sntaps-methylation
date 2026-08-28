#!/bin/bash
set -euo pipefail

BASE=/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation
cd "$BASE"

PY=/share/org/YZWL/yzwl_wujx/.conda/envs/snTAPS_analysis/bin/python
SCRIPT=03_scripts/34_umap_100kb_qcregress_fast.py
BIN_COUNTS=08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.bin_counts.tsv.gz
CORE_DIR=08_qc_summary/final_single_sample/cell_filter_strict_qc_core

LOGDIR=09_logs/region_umap_paper/qcregress_AB
CSUBDIR=03_scripts/region_matrix_paper/qcregress_AB_csub
mkdir -p "$LOGDIR" "$CSUBDIR"

echo "[INFO] PY=$PY"
echo "[INFO] SCRIPT=$SCRIPT"
echo "[INFO] BIN_COUNTS=$BIN_COUNTS"

if [ ! -x "$PY" ]; then
  echo "[ERROR] python not found: $PY" >&2
  exit 1
fi

if [ ! -s "$BIN_COUNTS" ]; then
  echo "[ERROR] bin counts not found: $BIN_COUNTS" >&2
  exit 1
fi

QCS=(
$(ls "$CORE_DIR"/*q20_q80*core5000*.global_per_cb.tsv 2>/dev/null | head -1 || true)
$(ls "$CORE_DIR"/*q20_q80*core2000*.global_per_cb.tsv 2>/dev/null | head -1 || true)
$(ls "$CORE_DIR"/*q25_q75*core5000*.global_per_cb.tsv 2>/dev/null | head -1 || true)
$(ls "$CORE_DIR"/*q25_q75*core2000*.global_per_cb.tsv 2>/dev/null | head -1 || true)
)

FEATURE_NAMES=(
all100kb
gene_body100kb
gene_body_no_promoter100kb
)

FEATURE_FILES=(
NONE
08_qc_summary/final_single_sample/feature_bins/gene_body_100kb_bin_ids.txt
08_qc_summary/final_single_sample/feature_bins/gene_body_no_promoter_100kb_bin_ids.txt
)

submit_one () {
  local qc="$1"
  local feat_name="$2"
  local feat_file="$3"
  local route="$4"
  local mode="$5"
  local qc_keys="$6"

  [ -s "$qc" ] || return 0

  local celltag
  celltag=$(basename "$qc" .global_per_cb.tsv)
  celltag=${celltag#snTAPS-FL-1.}

  local min_bin_cells=100
  if [[ "$celltag" == *core5000* ]]; then
    min_bin_cells=150
  elif [[ "$celltag" == *core2000* ]]; then
    min_bin_cells=60
  fi

  local tag="${feat_name}_${celltag}_${route}"
  local outdir="08_qc_summary/final_single_sample/region_umap_paper/100kb_qcregress_AB/${tag}"
  local prefix="snTAPS-FL-1.${tag}"

  local feature_arg=""
  if [ "$feat_file" != "NONE" ]; then
    if [ ! -s "$feat_file" ]; then
      echo "[WARN] missing feature file, skip: $feat_file"
      return 0
    fi
    feature_arg="--feature-bin-list $feat_file"
  fi

  local csub="$CSUBDIR/${tag}.csub"

  cat > "$csub" <<EOF
#!/bin/bash
#CSUB -J qAB_${feat_name}_${route}
#CSUB -q c02
#CSUB -n 16
#CSUB -R 'span[hosts=1]'
#CSUB -cwd $BASE
#CSUB -o $LOGDIR/${tag}.out
#CSUB -e $LOGDIR/${tag}.err

set -euo pipefail
cd $BASE

echo "[INFO] start \$(date)"
echo "[INFO] qc=$qc"
echo "[INFO] feature=$feat_name"
echo "[INFO] route=$route"
echo "[INFO] mode=$mode"
echo "[INFO] qc_keys=$qc_keys"
echo "[INFO] min_bin_cells=$min_bin_cells"

$PY $SCRIPT \\
  --global-qc "$qc" \\
  --bin-counts "$BIN_COUNTS" \\
  --outdir "$outdir" \\
  --prefix "$prefix" \\
  --mode "$mode" \\
  $feature_arg \\
  --cellbin-min-cov 1 \\
  --min-bin-cells $min_bin_cells \\
  --top-variable-bins 3000 \\
  --min-selected-bins 10 \\
  --n-pcs 30 \\
  --n-clusters 8 \\
  --n-neighbors 30 \\
  --min-dist 0.3 \\
  --n-components 2 \\
  --n-jobs 16 \\
  --seed 1 \\
  --qc-keys "$qc_keys"

echo "[INFO] done \$(date)"
EOF

  echo "[INFO] submit $tag"
  csub < "$csub"
}

for qc in "${QCS[@]}"; do
  [ -n "${qc:-}" ] || continue

  for i in "${!FEATURE_NAMES[@]}"; do
    feat_name="${FEATURE_NAMES[$i]}"
    feat_file="${FEATURE_FILES[$i]}"

    # 路线A：residual已经扣过全局modified_rate，所以只回归reads_used
    submit_one "$qc" "$feat_name" "$feat_file" "residual_depthOnly" "residual" "reads_used"

    # 路线B：seurat_norm不预先扣全局modified_rate，所以回归reads_used和modified_rate
    submit_one "$qc" "$feat_name" "$feat_file" "seurat_depthMod" "seurat_norm" "reads_used,modified_rate"
  done
done

echo "[INFO] submitted qcregress AB jobs."
