#!/usr/bin/env python3
import argparse
from collections import defaultdict
import gzip
import pysam
import pandas as pd

def load_regions(bed):
    regions_by_chrom = defaultdict(list)
    with open(bed) as f:
        for line in f:
            if not line.strip():
                continue
            a = line.rstrip("\n").split("\t")
            chrom = a[0]
            start = int(a[1])
            end = int(a[2])
            rid = a[3] if len(a) > 3 else f"{chrom}:{start}-{end}"
            regions_by_chrom[chrom].append((start, end, rid))
    for chrom in regions_by_chrom:
        regions_by_chrom[chrom].sort()
    return regions_by_chrom

def load_whitelist(path):
    s = set()
    with open(path) as f:
        for line in f:
            x = line.strip()
            if x:
                s.add(x.split()[0])
    return s

def get_cb(read):
    try:
        return read.get_tag("CB").split()[0]
    except Exception:
        return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam-list", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--whitelist", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--min-region-cov", type=int, default=5)
    args = ap.parse_args()

    whitelist = load_whitelist(args.whitelist)
    regions_by_chrom = load_regions(args.regions)
    fasta = pysam.FastaFile(args.fasta)

    with open(args.bam_list) as f:
        bams = [x.strip() for x in f if x.strip()]

    # counts[(CB, region_id)] = [modified, covered]
    counts = defaultdict(lambda: [0, 0])
    reads_used = 0

    for bam_path in bams:
        print(f"[INFO] BAM {bam_path}", flush=True)
        bam = pysam.AlignmentFile(bam_path, "rb")

        for chrom, regs in regions_by_chrom.items():
            if chrom not in bam.references:
                continue

            for start, end, rid in regs:
                ref = fasta.fetch(chrom, max(0, start - 2), end + 2).upper()
                offset = max(0, start - 2)

                for read in bam.fetch(chrom, start, end):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    if read.mapping_quality < args.min_mapq:
                        continue

                    cb = get_cb(read)
                    if cb not in whitelist:
                        continue

                    qseq = read.query_sequence
                    qqual = read.query_qualities
                    if qseq is None or qqual is None:
                        continue

                    reads_used += 1

                    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                        if qpos is None or rpos is None:
                            continue
                        if rpos < start or rpos >= end:
                            continue
                        if qqual[qpos] < args.min_baseq:
                            continue

                        idx = rpos - offset
                        if idx < 0 or idx >= len(ref):
                            continue

                        rb = ref[idx]
                        qb = qseq[qpos].upper()

                        # CpG C on forward strand
                        if rb == "C" and idx + 1 < len(ref) and ref[idx + 1] == "G":
                            counts[(cb, rid)][1] += 1
                            if qb == "T":
                                counts[(cb, rid)][0] += 1

                        # CpG C on reverse strand, represented as reference G preceded by C
                        elif rb == "G" and idx - 1 >= 0 and ref[idx - 1] == "C":
                            counts[(cb, rid)][1] += 1
                            if qb == "A":
                                counts[(cb, rid)][0] += 1

        bam.close()

    out_matrix = args.out_prefix + ".region_counts.tsv.gz"
    with gzip.open(out_matrix, "wt") as out:
        out.write("CB\tregion_id\tmodified_calls\tcovered_calls\tmodified_rate\n")
        for (cb, rid), (m, c) in counts.items():
            if c >= args.min_region_cov:
                out.write(f"{cb}\t{rid}\t{m}\t{c}\t{m/c:.6f}\n")

    summary = args.out_prefix + ".summary.tsv"
    with open(summary, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"bam_files\t{len(bams)}\n")
        out.write(f"whitelist_cells\t{len(whitelist)}\n")
        out.write(f"regions_total\t{sum(len(v) for v in regions_by_chrom.values())}\n")
        out.write(f"reads_used\t{reads_used}\n")
        out.write(f"cell_region_pairs\t{len(counts)}\n")

    print("[INFO] done", flush=True)
    print("[INFO] output", out_matrix, flush=True)

if __name__ == "__main__":
    main()
