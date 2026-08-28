#!/bin/bash
set -euo pipefail

PROJECT=/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation
RAWDIR=${PROJECT}/01_rawdata
LISTDIR=${PROJECT}/08_qc_summary/fastq_lists
OUT=${PROJECT}/08_qc_summary/sample_fastq_counts.tsv

mkdir -p "${LISTDIR}"

echo -e "sample\tR1_count\tR2_count" > "${OUT}"

find "${RAWDIR}" -mindepth 1 -maxdepth 1 -type d | sort | while read -r sample_dir
do
    sample=$(basename "$sample_dir")

    find "$sample_dir" -name "*_1.fq.gz" | sort > "${LISTDIR}/${sample}.R1.list"
    find "$sample_dir" -name "*_2.fq.gz" | sort > "${LISTDIR}/${sample}.R2.list"

    n1=$(wc -l < "${LISTDIR}/${sample}.R1.list")
    n2=$(wc -l < "${LISTDIR}/${sample}.R2.list")

    echo -e "${sample}\t${n1}\t${n2}" >> "${OUT}"
done

cat "${OUT}"
