#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: bash workflow/run_final_analysis.sh {preflight|promoter-matrix|final-umap|marker|figures} CONFIG" >&2
}

stage="${1:-}"
config="${2:-}"
if [[ -z "$stage" || -z "$config" ]]; then
  usage
  exit 2
fi
if [[ ! -f "$config" ]]; then
  echo "[ERROR] config not found: $config" >&2
  exit 2
fi

# shellcheck disable=SC1090
source "$config"

required=(PROJECT_ROOT SCRIPT_ROOT CONDA_INIT ANALYSIS_ENV REFERENCE_FASTA PROMOTER_BED BAM_LIST STRICT2_WHITELIST GLOBAL_QC PROMOTER_MATRIX_PREFIX P7_NO01_MATRIX FINAL_UMAP_OUTDIR FINAL_UMAP_PREFIX THREADS)
for name in "${required[@]}"; do
  if [[ -z "${!name:-}" ]]; then
    echo "[ERROR] missing config variable: $name" >&2
    exit 2
  fi
done

check_file() {
  [[ -s "$1" ]] || { echo "[ERROR] missing/non-empty file: $1" >&2; exit 1; }
}

preflight() {
  check_file "$CONDA_INIT"
  check_file "$REFERENCE_FASTA"
  check_file "$PROMOTER_BED"
  check_file "$BAM_LIST"
  check_file "$STRICT2_WHITELIST"
  check_file "$GLOBAL_QC"
  check_file "$SCRIPT_ROOT/region_matrix/40_call_region_matrix_by_cb.py"
  check_file "$SCRIPT_ROOT/region_matrix_paper/44_umap_region_matrix_paper_fast.py"
  check_file "$SCRIPT_ROOT/annotation/101_A_no01_rerun_marker_heatmap.py"
  check_file "$SCRIPT_ROOT/figures/110_make_A_main_figures.py"
  echo "[OK] stable final-analysis inputs and scripts are present"
  if [[ ! -s "$P7_NO01_MATRIX" ]]; then
    echo "[WARN] P7_no01 matrix is absent; recover/rebuild the missing matrix-generation step before final-umap" >&2
  fi
}

activate_analysis() {
  # shellcheck disable=SC1090
  source "$CONDA_INIT"
  conda activate "$ANALYSIS_ENV"
  cd "$PROJECT_ROOT"
}

case "$stage" in
  preflight)
    preflight
    ;;
  promoter-matrix)
    preflight
    activate_analysis
    mkdir -p "$(dirname "$PROMOTER_MATRIX_PREFIX")"
    python3 "$SCRIPT_ROOT/region_matrix/40_call_region_matrix_by_cb.py" \
      --bam-list "$BAM_LIST" \
      --regions "$PROMOTER_BED" \
      --fasta "$REFERENCE_FASTA" \
      --whitelist "$STRICT2_WHITELIST" \
      --out-prefix "$PROMOTER_MATRIX_PREFIX" \
      --min-mapq 30 \
      --min-baseq 20 \
      --min-region-cov 5
    ;;
  final-umap)
    preflight
    check_file "$P7_NO01_MATRIX"
    activate_analysis
    mkdir -p "$FINAL_UMAP_OUTDIR"
    python3 "$SCRIPT_ROOT/region_matrix_paper/44_umap_region_matrix_paper_fast.py" \
      --region-counts "$P7_NO01_MATRIX" \
      --global-qc "$GLOBAL_QC" \
      --outdir "$FINAL_UMAP_OUTDIR" \
      --prefix "$FINAL_UMAP_PREFIX" \
      --mode residual \
      --max-cells 0 \
      --min-cell-regions 500 \
      --min-region-cells 1000 \
      --top-variable-regions 2000 \
      --regress-qc \
      --drop-pc1 \
      --n-pcs 10 \
      --n-neighbors 100 \
      --min-dist 0.4 \
      --n-clusters 6 \
      --n-jobs "$THREADS" \
      --seed 1
    ;;
  marker)
    preflight
    activate_analysis
    python3 "$SCRIPT_ROOT/annotation/101_A_no01_rerun_marker_heatmap.py"
    ;;
  figures)
    preflight
    activate_analysis
    python3 "$SCRIPT_ROOT/figures/110_make_A_main_figures.py"
    ;;
  *)
    usage
    exit 2
    ;;
esac
