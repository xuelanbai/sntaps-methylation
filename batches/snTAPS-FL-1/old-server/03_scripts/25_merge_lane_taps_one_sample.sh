#!/bin/bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: bash 03_scripts/25_merge_lane_taps_one_sample.sh SAMPLE"
    exit 1
fi

SAMPLE="$1"
PROJECT=/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation
OUTDIR=${PROJECT}/08_qc_summary/downstream_test/${SAMPLE}
LANEDIR=${OUTDIR}/by_lane
MERGEDIR=${OUTDIR}/merged_full
TMPDIR=${OUTDIR}/tmp_merge_${SAMPLE}

mkdir -p "$MERGEDIR" "$TMPDIR"

cd "$PROJECT"

echo "[INFO] sample: $SAMPLE"
echo "[INFO] lanedir: $LANEDIR"

n_summary=$(find "$LANEDIR" -name "*.full.summary.tsv" | wc -l)
n_global=$(find "$LANEDIR" -name "*.full.global_per_cb.tsv" | wc -l)
n_bins=$(find "$LANEDIR" -name "*.full.bin_counts.tsv.gz" | wc -l)

echo "[INFO] summary files: $n_summary"
echo "[INFO] global files: $n_global"
echo "[INFO] bin files: $n_bins"

if [ "$n_summary" -ne 8 ] || [ "$n_global" -ne 8 ] || [ "$n_bins" -ne 8 ]; then
    echo "[ERROR] expected 8 lane outputs for $SAMPLE"
    exit 1
fi

# 1) merge bin_counts
echo "[INFO] merging bin_counts: $(date)"

(
    echo -e "CB\tbin_id\tmodified_calls\tcovered_calls\tmodified_rate"

    zcat ${LANEDIR}/${SAMPLE}_L*.full.bin_counts.tsv.gz \
      | awk 'BEGIN{OFS="\t"} $1=="CB" && $2=="bin_id"{next} {print $1,$2,$3,$4}' \
      | sort -T "$TMPDIR" -S 20G --parallel=8 -k1,1 -k2,2 \
      | awk '
        BEGIN{FS=OFS="\t"}
        {
            key=$1"\t"$2
            if (NR==1) {
                prev=key
                cb=$1
                bin=$2
                mod=$3
                cov=$4
            } else if (key==prev) {
                mod+=$3
                cov+=$4
            } else {
                rate=(cov>0 ? mod/cov : 0)
                print cb,bin,mod,cov,sprintf("%.6f",rate)
                prev=key
                cb=$1
                bin=$2
                mod=$3
                cov=$4
            }
        }
        END{
            if (NR>0) {
                rate=(cov>0 ? mod/cov : 0)
                print cb,bin,mod,cov,sprintf("%.6f",rate)
            }
        }'
) | gzip -1 > ${MERGEDIR}/${SAMPLE}.full_merged.bin_counts.tsv.gz

echo "[INFO] wrote ${MERGEDIR}/${SAMPLE}.full_merged.bin_counts.tsv.gz"

# 2) sum reads_used per CB from lane global files
echo "[INFO] merging reads_used per CB"

awk '
BEGIN{FS=OFS="\t"}
FNR==1{next}
{
    reads[$1]+=$2
}
END{
    print "CB","reads_used"
    for (cb in reads) print cb,reads[cb]
}' ${LANEDIR}/${SAMPLE}_L*.full.global_per_cb.tsv \
  | sort -k1,1 \
  > ${TMPDIR}/${SAMPLE}.reads_by_cb.tsv

# 3) calculate global_per_cb from merged bin_counts
echo "[INFO] building global_per_cb"

zcat ${MERGEDIR}/${SAMPLE}.full_merged.bin_counts.tsv.gz \
  | awk '
    BEGIN{FS=OFS="\t"}
    NR==1{next}
    {
        cb=$1
        mod[cb]+=$3
        cov[cb]+=$4
        bins[cb]++
    }
    END{
        print "CB","modified_calls","covered_calls","modified_rate","covered_bins"
        for (cb in cov) {
            rate=(cov[cb]>0 ? mod[cb]/cov[cb] : 0)
            print cb,mod[cb],cov[cb],sprintf("%.6f",rate),bins[cb]
        }
    }' \
  | sort -k1,1 \
  > ${TMPDIR}/${SAMPLE}.global_cov_by_cb.tsv

# 4) join reads_used and global cov
join -t $'\t' -a 1 -e 0 -o 1.1,1.2,2.2,2.3,2.4,2.5 \
  <(tail -n +2 ${TMPDIR}/${SAMPLE}.reads_by_cb.tsv | sort -k1,1) \
  <(tail -n +2 ${TMPDIR}/${SAMPLE}.global_cov_by_cb.tsv | sort -k1,1) \
  | awk 'BEGIN{FS=OFS="\t"; print "CB","reads_used","modified_calls","covered_calls","modified_rate","covered_bins"} {print}' \
  > ${MERGEDIR}/${SAMPLE}.full_merged.global_per_cb.tsv

# 5) summary
echo "[INFO] building summary"

reads_used=$(awk 'FNR>1{sum+=$2} END{print sum+0}' ${LANEDIR}/${SAMPLE}_L*.full.global_per_cb.tsv)
cells=$(awk 'NR>1{n++} END{print n+0}' ${MERGEDIR}/${SAMPLE}.full_merged.global_per_cb.tsv)
total_mod=$(awk 'NR>1{sum+=$3} END{print sum+0}' ${MERGEDIR}/${SAMPLE}.full_merged.global_per_cb.tsv)
total_cov=$(awk 'NR>1{sum+=$4} END{print sum+0}' ${MERGEDIR}/${SAMPLE}.full_merged.global_per_cb.tsv)
rate=$(awk -v m="$total_mod" -v c="$total_cov" 'BEGIN{printf "%.6f", (c>0?m/c:0)}')

{
    echo -e "metric\tvalue"
    echo -e "sample\t${SAMPLE}"
    echo -e "lane_files\t8"
    echo -e "reads_used\t${reads_used}"
    echo -e "cells_with_covered_CpG\t${cells}"
    echo -e "total_modified_calls\t${total_mod}"
    echo -e "total_covered_calls\t${total_cov}"
    echo -e "overall_modified_rate\t${rate}"
} > ${MERGEDIR}/${SAMPLE}.full_merged.summary.tsv

echo "[INFO] wrote ${MERGEDIR}/${SAMPLE}.full_merged.global_per_cb.tsv"
echo "[INFO] wrote ${MERGEDIR}/${SAMPLE}.full_merged.summary.tsv"

rm -rf "$TMPDIR"

echo "[INFO] done: $(date)"
