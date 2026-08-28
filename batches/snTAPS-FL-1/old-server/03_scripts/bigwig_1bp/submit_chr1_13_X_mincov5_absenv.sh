#!/bin/bash
set -eo pipefail

BASE=/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation
PY=/share/org/YZWL/yzwl_wujx/.conda/envs/snTAPS_analysis/bin/python
BWTOOL=/share/org/YZWL/yzwl_wujx/.conda/envs/foodie_env/bin/bedGraphToBigWig

SCRIPT=${BASE}/03_scripts/bigwig_1bp/convert_one_chr_cpg_to_bw.py
SITE_DIR=${BASE}/08_qc_summary/final_single_sample/site_conversion
OUTROOT=${BASE}/08_qc_summary/final_single_sample/bigwig_1bp_methyl
GENOME=${OUTROOT}/genome.sizes
LOGDIR=${BASE}/09_logs/bigwig_1bp

mkdir -p "$OUTROOT" "$LOGDIR"

for CHR in 1 13 X
do
  IN=${SITE_DIR}/snTAPS-FL-1.${CHR}.CpG_site_conversion.tsv.gz
  OUTDIR=${OUTROOT}/bw_chr_${CHR}_mincov5
  CSUB=${LOGDIR}/bw_chr_${CHR}.csub

  cat > "$CSUB" <<EOF
#!/bin/bash
#CSUB -J bw_${CHR}_m5
#CSUB -q c02
#CSUB -n 4
#CSUB -R 'span[hosts=1]'
#CSUB -cwd ${BASE}
#CSUB -o ${LOGDIR}/bw_chr_${CHR}.out
#CSUB -e ${LOGDIR}/bw_chr_${CHR}.err

set -eo pipefail

cd ${BASE}

PY=${PY}
BWTOOL=${BWTOOL}
SCRIPT=${SCRIPT}
IN=${IN}
OUTDIR=${OUTDIR}
GENOME=${GENOME}

echo "[INFO] host=\$(hostname)"
echo "[INFO] python=\$PY"
echo "[INFO] bwtool=\$BWTOOL"
echo "[INFO] input=\$IN"
echo "[INFO] outdir=\$OUTDIR"

test -s "\$PY"
test -s "\$BWTOOL"
test -s "\$SCRIPT"
test -s "\$IN"
test -s "\$GENOME"

"\$PY" - <<'PYCHK'
import pandas as pd
print("[INFO] pandas", pd.__version__)
PYCHK

mkdir -p "\$OUTDIR"

"\$PY" "\$SCRIPT" \
  --input "\$IN" \
  --chrom "${CHR}" \
  --outdir "\$OUTDIR" \
  --genome-size "\$GENOME" \
  --mincov 5 \
  --chunk-size 10000000 \
  --bedGraphToBigWig "\$BWTOOL"

echo "[DONE] chr${CHR}"
find "\$OUTDIR" -maxdepth 1 -type f -printf "%p\t%s bytes\n" | sort
EOF

  echo "[SUBMIT] $CSUB"
  csub < "$CSUB"
done
