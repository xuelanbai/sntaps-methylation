#!/bin/bash
set -euo pipefail

PROJECT=/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation
RAWROOT=${PROJECT}/00_raw_release/zhaoqiulin0430017/X101SC26014480-Z01-J019_multipath
OUTDIR=${PROJECT}/01_rawdata
SUMMARY=${PROJECT}/08_qc_summary/raw_fastq_link_summary.tsv

mkdir -p "${OUTDIR}" "${PROJECT}/08_qc_summary"

echo -e "package\tsample\tlane\tread\tfastq" > "${SUMMARY}"

find "${RAWROOT}" -path "*/01.RawData/*/*.fq.gz" \
  ! -path "*/Undetermined/*" \
  | sort \
  | while read -r fq
do
    base=$(basename "$fq")
    package=$(echo "$fq" | grep -oE 'X101SC26014480-Z01-J019_[0-9]+' | head -1)
    sample=$(basename "$(dirname "$fq")")
    lane=$(echo "$base" | grep -oE '_L[0-9]+_' | tr -d '_')
    readnum=$(echo "$base" | sed -E 's/.*_([12])\.fq\.gz$/\1/')
    read="R${readnum}"

    mkdir -p "${OUTDIR}/${sample}"
    ln -sfn "$fq" "${OUTDIR}/${sample}/${base}"

    echo -e "${package}\t${sample}\t${lane}\t${read}\t${OUTDIR}/${sample}/${base}" >> "${SUMMARY}"
done

echo "[INFO] linked fastq number:"
find "${OUTDIR}" -maxdepth 2 -name "*.fq.gz" | wc -l

echo "[INFO] sample list:"
cut -f2 "${SUMMARY}" | tail -n +2 | sort -u
