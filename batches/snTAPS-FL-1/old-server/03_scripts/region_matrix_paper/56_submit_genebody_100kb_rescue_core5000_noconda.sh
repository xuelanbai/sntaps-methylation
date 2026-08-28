#!/bin/bash
set -euo pipefail

BASE=/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation
cd "$BASE"

PYTHON_BIN="/share/org/YZWL/yzwl_wujx/.conda/envs/snTAPS_analysis/bin/python3"

QC_DIR="08_qc_summary/final_single_sample/cell_filter_strict_qc_core"

GLOBAL_QC=$(ls ${QC_DIR}/*q25_q75_ge100k*core5000*global_per_cb.tsv 2>/dev/null | head -n 1 || true)
if [ -z "${GLOBAL_QC}" ]; then
  GLOBAL_QC=$(ls ${QC_DIR}/*q20_q80_ge100k*core5000*global_per_cb.tsv 2>/dev/null | head -n 1 || true)
fi

if [ -z "${GLOBAL_QC}" ] || [ ! -s "${GLOBAL_QC}" ]; then
  echo "[ERROR] cannot find strict core5000 global_per_cb.tsv in ${QC_DIR}"
  ls -lh ${QC_DIR}/*global_per_cb.tsv || true
  exit 1
fi

BIN_COUNTS="08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.bin_counts.tsv.gz"

declare -A FEATURE_LISTS
FEATURE_LISTS[gene_body]="08_qc_summary/final_single_sample/feature_bins/gene_body_100kb_bin_ids.txt"
FEATURE_LISTS[gene_body_no_promoter]="08_qc_summary/final_single_sample/feature_bins/gene_body_no_promoter_100kb_bin_ids.txt"

MODES=("raw" "residual" "seurat_norm")

LOGDIR="09_logs/region_umap_paper/100kb_rescue_genebody_core5000"
mkdir -p "$LOGDIR"

for FEATURE in gene_body gene_body_no_promoter; do
  FEATURE_BIN_LIST="${FEATURE_LISTS[$FEATURE]}"

  if [ ! -s "$FEATURE_BIN_LIST" ]; then
    echo "[ERROR] missing feature bin list: $FEATURE_BIN_LIST"
    exit 1
  fi

  for MODE in "${MODES[@]}"; do
    OUTDIR="08_qc_summary/final_single_sample/region_umap_paper/100kb_rescue/${FEATURE}100kb_${MODE}_strict_core5000_cellbinCov1_top3000"
    PREFIX="snTAPS-FL-1.${FEATURE}100kb_${MODE}_strict_core5000_cellbinCov1_top3000"
    CSUB_SCRIPT="${LOGDIR}/${PREFIX}.csub"

    cat > "$CSUB_SCRIPT" <<EOS
#!/bin/bash
#CSUB -J gb100kb_${FEATURE}_${MODE}
#CSUB -q c02
#CSUB -n 16
#CSUB -R 'span[hosts=1]'
#CSUB -cwd $BASE
#CSUB -o ${LOGDIR}/${PREFIX}.out
#CSUB -e ${LOGDIR}/${PREFIX}.err

set -euo pipefail
cd $BASE

echo "[INFO] PYTHON_BIN=$PYTHON_BIN"
$PYTHON_BIN -V

echo "[INFO] GLOBAL_QC=${GLOBAL_QC}"
echo "[INFO] BIN_COUNTS=${BIN_COUNTS}"
echo "[INFO] FEATURE=${FEATURE}"
echo "[INFO] FEATURE_BIN_LIST=${FEATURE_BIN_LIST}"
echo "[INFO] MODE=${MODE}"
echo "[INFO] OUTDIR=${OUTDIR}"
echo "[INFO] PREFIX=${PREFIX}"

$PYTHON_BIN 03_scripts/33_umap_100kb_modes.py \
  --global-qc "${GLOBAL_QC}" \
  --bin-counts "${BIN_COUNTS}" \
  --outdir "${OUTDIR}" \
  --prefix "${PREFIX}" \
  --mode "${MODE}" \
  --feature-bin-list "${FEATURE_BIN_LIST}" \
  --min-reads 0 \
  --min-covered-calls 0 \
  --min-covered-bins 0 \
  --cellbin-min-cov 1 \
  --max-cells 200000 \
  --min-bin-cells 300 \
  --top-variable-bins 3000 \
  --n-pcs 30 \
  --n-clusters 8 \
  --n-jobs 16 \
  --seed 1
EOS

    echo "[INFO] submit $CSUB_SCRIPT"
    csub < "$CSUB_SCRIPT"
  done
done

echo "[INFO] submitted genebody 100kb rescue jobs."
