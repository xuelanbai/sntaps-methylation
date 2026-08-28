#!/usr/bin/env python3
import argparse
from collections import defaultdict
import gzip
import pysam

def get_cb(read):
    try:
        return read.get_tag("CB").split()[0]
    except Exception:
        return None

def load_whitelist(path):
    if not path:
        return None
    s = set()
    with open(path) as f:
        for line in f:
            x = line.strip()
            if x and x != "CB":
                s.add(x.split()[0])
    return s

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam-list", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--whitelist", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--min-baseq", type=int, default=20)
    args = ap.parse_args()

    whitelist = load_whitelist(args.whitelist)

    fasta = pysam.FastaFile(args.fasta)
    ref = fasta.fetch(args.chrom).upper()
    chrom_len = len(ref)

    bams = []
    with open(args.bam_list) as f:
        for line in f:
            line = line.strip()
            if line:
                bams.append(line)

    # key = CpG C position, 0-based
    cov = defaultdict(int)
    mod = defaultdict(int)

    total_reads = 0
    used_reads = 0

    for bam_path in bams:
        print(f"[INFO] chrom={args.chrom} bam={bam_path}", flush=True)
        bam = pysam.AlignmentFile(bam_path, "rb")

        for read in bam.fetch(args.chrom):
            total_reads += 1

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < args.min_mapq:
                continue

            if whitelist is not None:
                cb = get_cb(read)
                if cb not in whitelist:
                    continue

            qseq = read.query_sequence
            qqual = read.query_qualities
            if qseq is None or qqual is None:
                continue

            used_reads += 1

            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if qpos is None or rpos is None:
                    continue
                if qqual[qpos] < args.min_baseq:
                    continue
                if rpos < 0 or rpos >= chrom_len:
                    continue

                rb = ref[rpos]
                qb = qseq[qpos].upper()

                # Forward strand CpG C: C followed by G
                if rb == "C" and rpos + 1 < chrom_len and ref[rpos + 1] == "G":
                    site = rpos
                    cov[site] += 1
                    if qb == "T":
                        mod[site] += 1

                # Reverse strand CpG C is reference G preceded by C.
                # Map this observation to the same CpG C coordinate: rpos - 1
                elif rb == "G" and rpos - 1 >= 0 and ref[rpos - 1] == "C":
                    site = rpos - 1
                    cov[site] += 1
                    if qb == "A":
                        mod[site] += 1

        bam.close()

    with gzip.open(args.out, "wt") as out:
        out.write("chrom\tstart\tend\tmodified_calls\tcovered_calls\tmodified_rate\n")
        for site in sorted(cov):
            c = cov[site]
            m = mod[site]
            rate = m / c if c else 0
            # CpG interval covers two bases: C and G
            out.write(f"{args.chrom}\t{site}\t{site+2}\t{m}\t{c}\t{rate:.6f}\n")

    print(f"[INFO] done chrom={args.chrom} sites={len(cov)} total_reads={total_reads} used_reads={used_reads}", flush=True)

if __name__ == "__main__":
    main()
