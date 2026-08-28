#!/bin/bash
set -eo pipefail

cd /share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation

SITE_DIR=08_qc_summary/final_single_sample/site_conversion
OUTROOT=08_qc_summary/final_single_sample/bigwig_1bp_methyl
GENOME_SIZE=${OUTROOT}/genome.sizes
LOGDIR=09_logs/bigwig_1bp

mkdir -p "$LOGDIR"

for chr in 1 13 X
do
    input=${SITE_DIR}/snTAPS-FL-1.${chr}.CpG_site_conversion.tsv.gz
    outdir=${OUTROOT}/bw_chr_${chr}_mincov5
    csub_file=${LOGDIR}/bw_chr_${chr}.csub

    cat > "$csub_file" <<EOF
#!/bin/bash
#CSUB -J bw_chr_${chr}
#CSUB -q c02
#CSUB -n 8
#CSUB -R 'span[hosts=1]'
#CSUB -cwd /share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation
#CSUB -o ${LOGDIR}/bw_chr_${chr}.out
#CSUB -e ${LOGDIR}/bw_chr_${chr}.err

set -eo pipefail

cd /share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation

source ~/.bashrc || true
conda activate snTAPS_analysis || true

set -u

BWTOOL=\$(command -v bedGraphToBigWig || true)
if [[ -z "\$BWTOOL" ]]; then
    BWTOOL=\$(find /share /opt /usr/local -name bedGraphToBigWig 2>/dev/null | head -1 || true)
fi
if [[ -z "\$BWTOOL" ]]; then
    echo "[ERROR] bedGraphToBigWig not found"
    exit 1
fi

mkdir -p "$outdir"

python3 03_scripts/bigwig_1bp/convert_one_chr_cpg_to_bw.py \\
  --input "$input" \\
  --genome-size "$GENOME_SIZE" \\
  --outdir "$outdir" \\
  --min-cov 5 \\
  --chunk-size 10000000 \\
  --bedGraphToBigWig "\$BWTOOL"

echo "[DONE] chr${chr}"
EOF

    csub < "$csub_file"
done
