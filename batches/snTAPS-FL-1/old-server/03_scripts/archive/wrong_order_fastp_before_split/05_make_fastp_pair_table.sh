#!/bin/bash
set -euo pipefail

PROJECT=/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation
RAWDIR=${PROJECT}/01_rawdata
OUT=${PROJECT}/08_qc_summary/fastp_pair_table.tsv

mkdir -p ${PROJECT}/08_qc_summary

echo -e "sample\tlane\tR1\tR2" > "${OUT}"

find "${RAWDIR}" -mindepth 1 -maxdepth 1 -type d | sort | while read -r sample_dir
do
    sample=$(basename "$sample_dir")

    find "$sample_dir" -name "*_1.fq.gz" | sort | while read -r r1
    do
        r2=$(echo "$r1" | sed 's/_1\.fq\.gz$/_2.fq.gz/')
        lane=$(basename "$r1" | grep -oE '_L[0-9]+_' | tr -d '_')

        if [ ! -e "$r2" ]; then
            echo "[ERROR] missing R2 for $r1" >&2
            exit 1
        fi

        echo -e "${sample}\t${lane}\t${r1}\t${r2}" >> "${OUT}"
    done
done

cat "${OUT}"
