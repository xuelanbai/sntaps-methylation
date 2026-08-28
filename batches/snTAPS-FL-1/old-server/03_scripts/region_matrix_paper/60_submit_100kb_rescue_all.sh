#!/bin/bash
set -euo pipefail

BASE=/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation
cd "$BASE"

GLOBAL_QC="08_qc_summary/final_single_sample/cell_filter_strict_q10_q95_ge100k/snTAPS-FL-1.strict_q10_q95_ge100k.global_per_cb.tsv"
SCRIPT="03_scripts/33_umap_100kb_modes.py"

# 使用最终合并后的per-cell 100kb bin_counts，不自动选择archive/downstream_test中间文件
BIN_COUNTS="08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.bin_counts.tsv.gz"

echo "[INFO] using BIN_COUNTS=$BIN_COUNTS"
echo "[INFO] using GLOBAL_QC=$GLOBAL_QC"

mkdir -p 09_logs/umap_100kb_rescue
mkdir -p 08_qc_summary/final_single_sample/region_umap_paper/100kb_rescue

run_one () {
  name=$1
  mode=$2
  feature=$3
  cov=$4
  topn=$5

  outdir="08_qc_summary/final_single_sample/region_umap_paper/100kb_rescue/${name}_${mode}_cellbinCov${cov}_top${topn}"
  prefix="snTAPS-FL-1.${name}_${mode}_q10_q95_ge100k_cellbinCov${cov}_top${topn}"

  job="umap_${name}_${mode}_c${cov}_t${topn}"
  shfile="09_logs/umap_100kb_rescue/${job}.csub"

  mkdir -p "$outdir"

  cat > "$shfile" <<EOS
#!/bin/bash
#CSUB -J ${job}
#CSUB -q c02
#CSUB -n 16
#CSUB -o 09_logs/umap_100kb_rescue/${job}.out
#CSUB -e 09_logs/umap_100kb_rescue/${job}.err

set -euo pipefail
cd "$BASE"

source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate snTAPS_analysis

cmd=(python3 "$SCRIPT"
  --global-qc "$GLOBAL_QC"
  --bin-counts "$BIN_COUNTS"
  --outdir "$outdir"
  --prefix "$prefix"
  --mode "$mode"
  --min-reads 0
  --min-covered-calls 0
  --min-covered-bins 0
  --cellbin-min-cov "$cov"
  --max-cells 200000
  --min-bin-cells 300
  --top-variable-bins "$topn"
  --n-pcs 30
  --n-clusters 8
  --n-jobs 16
  --seed 1
)

if [ "$feature" != "NONE" ]; then
  cmd+=(--feature-bin-list "$feature")
fi

echo "[INFO] \${cmd[@]}"
"\${cmd[@]}"
EOS

  csub < "$shfile"
}

# all 100kb bins：不限制feature
for mode in residual raw seurat_norm; do
  for cov in 1 5; do
    run_one "all100kb" "$mode" "NONE" "$cov" 3000
  done
done

# promoter / gene body / gene body no promoter / gene region 100kb bins
for feature_name in promoter gene_body gene_body_no_promoter gene_region; do
  feature_file="08_qc_summary/final_single_sample/feature_bins/${feature_name}_100kb_bin_ids.txt"
  if [ -s "$feature_file" ]; then
    for mode in residual raw seurat_norm; do
      for cov in 1 5; do
        run_one "${feature_name}100kb" "$mode" "$feature_file" "$cov" 3000
      done
    done
  else
    echo "[WARN] missing feature file: $feature_file"
  fi
done
