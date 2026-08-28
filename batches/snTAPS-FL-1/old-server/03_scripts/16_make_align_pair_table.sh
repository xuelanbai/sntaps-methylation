#!/bin/bash
set -euo pipefail

PROJECT=/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation
FQDIR=${PROJECT}/05_fastp/after_split
OUT=${PROJECT}/08_qc_summary/align_pair_table.tsv

echo -e "sample\tlane\tR1\tR2" > "$OUT"

find "$FQDIR" -mindepth 1 -maxdepth 1 -type d | sort | while read -r sample_dir
do
    sample=$(basename "$sample_dir")

    find "$sample_dir" -name "*.fastp_R1.fq.gz" | sort | while read -r r1
    do
        r2="${r1%.fastp_R1.fq.gz}.fastp_R2.fq.gz"
        lane=$(basename "$r1" | grep -oE '_L[0-9]+' | tr -d '_' || true)

        if [ ! -e "$r2" ]; then
            echo "[ERROR] missing R2: $r2" >&2
            exit 1
        fi

        echo -e "${sample}\t${lane}\t${r1}\t${r2}" >> "$OUT"
    done
done

cat "$OUT"
