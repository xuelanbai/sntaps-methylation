#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path

def open_text(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def load_chrom_sizes(path):
    chroms = []
    sizes = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            chrom, size = line.strip().split()[:2]
            size = int(size)
            chroms.append(chrom)
            sizes[chrom] = size
    return chroms, sizes

def make_bins(chroms, sizes, bin_size):
    bins = {}
    for chrom in chroms:
        nbin = (sizes[chrom] + bin_size - 1) // bin_size
        bins[chrom] = [0] * nbin
    return bins

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--bin-size", type=int, default=100000)
    ap.add_argument("--min-cpg", type=int, default=20)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    chroms, sizes = load_chrom_sizes(args.chrom_sizes)
    keep = set(chroms)
    cpg_counts = make_bins(chroms, sizes, args.bin_size)

    bins_bed = outdir / "pig_100kb_bins.bed"
    cpg_bed_gz = outdir / "pig_CpG_sites.bed.gz"
    cpg_count_bed = outdir / "pig_100kb_bins.cpg_count.bed"
    cpg20_bed = outdir / f"pig_100kb_bins.cpg{args.min_cpg}.bed"
    summary = outdir / "pig_cpg_bin_summary.txt"

    # 先输出100-kb bins
    total_bins = 0
    with open(bins_bed, "w") as fo:
        for chrom in chroms:
            size = sizes[chrom]
            for start in range(0, size, args.bin_size):
                end = min(start + args.bin_size, size)
                center = (start + end) // 2
                bin_id = f"{chrom}:{start}-{end}"
                fo.write(f"{chrom}\t{start}\t{end}\t{bin_id}\t0\t{center}\n")
                total_bins += 1

    # 扫描CpG并计数
    total_cpg = 0
    active = False
    chrom = None
    pos = 0
    prev_base = ""

    with open_text(args.fasta, "rt") as fi, gzip.open(cpg_bed_gz, "wt", compresslevel=1) as cpg_out:
        for line in fi:
            line = line.rstrip()
            if not line:
                continue

            if line.startswith(">"):
                chrom = line[1:].split()[0]
                active = chrom in keep
                pos = 0
                prev_base = ""
                continue

            seq = line.upper()

            if not active:
                pos += len(seq)
                prev_base = seq[-1] if seq else ""
                continue

            for i, base in enumerate(seq):
                if prev_base == "C" and base == "G":
                    cpg_start = pos + i - 1
                    cpg_end = pos + i + 1

                    if cpg_start >= 0:
                        cpg_out.write(f"{chrom}\t{cpg_start}\t{cpg_end}\n")
                        total_cpg += 1

                        b = cpg_start // args.bin_size
                        if b < len(cpg_counts[chrom]):
                            cpg_counts[chrom][b] += 1

                prev_base = base

            pos += len(seq)

    # 输出每个bin的CpG数量
    kept_bins = 0
    with open(cpg_count_bed, "w") as all_out, open(cpg20_bed, "w") as keep_out:
        for chrom in chroms:
            size = sizes[chrom]
            for idx, cpg_n in enumerate(cpg_counts[chrom]):
                start = idx * args.bin_size
                end = min(start + args.bin_size, size)
                center = (start + end) // 2
                bin_id = f"{chrom}:{start}-{end}"

                row = f"{chrom}\t{start}\t{end}\t{bin_id}\t{cpg_n}\t{center}\n"
                all_out.write(row)

                if cpg_n >= args.min_cpg:
                    keep_out.write(row)
                    kept_bins += 1

    with open(summary, "w") as fo:
        fo.write(f"fasta\t{args.fasta}\n")
        fo.write(f"chrom_sizes\t{args.chrom_sizes}\n")
        fo.write(f"bin_size\t{args.bin_size}\n")
        fo.write(f"min_cpg\t{args.min_cpg}\n")
        fo.write(f"primary_chromosomes\t{len(chroms)}\n")
        fo.write(f"total_bins\t{total_bins}\n")
        fo.write(f"total_cpg_sites\t{total_cpg}\n")
        fo.write(f"kept_bins_cpg_ge_{args.min_cpg}\t{kept_bins}\n")

    print(f"[INFO] wrote {bins_bed}")
    print(f"[INFO] wrote {cpg_bed_gz}")
    print(f"[INFO] wrote {cpg_count_bed}")
    print(f"[INFO] wrote {cpg20_bed}")
    print(f"[INFO] wrote {summary}")

if __name__ == "__main__":
    main()
