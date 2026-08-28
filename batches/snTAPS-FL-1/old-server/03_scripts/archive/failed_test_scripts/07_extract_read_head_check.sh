#!/bin/bash
set -euo pipefail

PROJECT=/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation
RAWDIR=${PROJECT}/01_rawdata
OUTDIR=${PROJECT}/08_qc_summary/read_structure

mkdir -p "${OUTDIR}"

R1_LIST=${OUTDIR}/test_R1.list

find "${RAWDIR}" -maxdepth 2 -name "*_1.fq.gz" | sort | head -8 > "${R1_LIST}"

echo "[INFO] R1 files selected:"
cat "${R1_LIST}"
echo "[INFO] R1 count:"
wc -l "${R1_LIST}"

if [ ! -s "${R1_LIST}" ]; then
    echo "[ERROR] no R1 fastq found in ${RAWDIR}" >&2
    exit 1
fi

while read -r r1
do
    r2="${r1%_1.fq.gz}_2.fq.gz"

    if [ ! -e "$r2" ]; then
        echo "[ERROR] missing R2: $r2" >&2
        exit 1
    fi

    sample=$(basename "$(dirname "$r1")")
    lane=$(basename "$r1" | grep -oE '_L[0-9]+_' | tr -d '_' || true)

    if [ -z "$lane" ]; then
        lane="unknownLane"
    fi

    out="${OUTDIR}/${sample}_${lane}.head_check.txt"
    tmp1="${OUTDIR}/${sample}_${lane}.R1.tmp"
    tmp2="${OUTDIR}/${sample}_${lane}.R2.tmp"

    zcat "$r1" | head -40 > "$tmp1"
    zcat "$r2" | head -40 > "$tmp2"

    {
        echo "### ${sample} ${lane}"
        echo "R1=${r1}"
        echo "R2=${r2}"
        echo

        paste "$tmp1" "$tmp2" | awk '
            NR%4==1 {
                print "HEADER_R1\t"$1
                print "HEADER_R2\t"$2
            }
            NR%4==2 {
                print "R1_first40\t"substr($1,1,40)
                print "R2_first40\t"substr($2,1,40)
                print ""
            }
        '
    } > "$out"

    rm -f "$tmp1" "$tmp2"

    echo "[INFO] wrote ${out}"
done < "${R1_LIST}"
