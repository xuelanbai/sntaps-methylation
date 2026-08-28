#!/usr/bin/env python3
import argparse
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path
import pysam

def load_whitelist(path):
    if not path:
        return None
    s = set()
    with open(path) as f:
        for line in f:
            x = line.strip().split()
            if x:
                s.add(x[0])
    return s

def load_bins(path):
    bin_map = defaultdict(dict)
    bin_size_guess = None
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            chrom = p[0]
            start = int(p[1])
            end = int(p[2])
            bin_id = p[3]
            idx = start // 100000
            bin_map[chrom][idx] = bin_id
            if bin_size_guess is None:
                bin_size_guess = end - start
    return bin_map

def get_cb(read):
    try:
        cb = read.get_tag("CB")
    except KeyError:
        return None

    # Some BAMs were generated with FASTQ comments like:
    # CB:Z:<barcode> CN:Z:<barcode_name>
    # pysam may read the whole string as CB value.
    # Keep only the pure barcode sequence before whitespace.
    cb = cb.split()[0]
    return cb

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam-list", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--bins", required=True)
    ap.add_argument("--whitelist", default=None)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--max-reads", type=int, default=0,
                    help="0 means no limit; limit is total reads passing read-level filters.")
    args = ap.parse_args()

    whitelist = load_whitelist(args.whitelist)
    bin_map = load_bins(args.bins)

    fa = pysam.FastaFile(args.fasta)

    cb_mod = Counter()
    cb_cov = Counter()
    cb_reads = Counter()
    cb_bin_mod = Counter()
    cb_bin_cov = Counter()

    summary = Counter()
    reads_used = 0

    with open(args.bam_list) as f:
        bam_paths = [x.strip() for x in f if x.strip()]

    for bam_path in bam_paths:
        sys.stderr.write(f"[INFO] processing {bam_path}\n")
        bam = pysam.AlignmentFile(bam_path, "rb")

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                summary["unmapped"] += 1
                continue
            if read.is_secondary or read.is_supplementary:
                summary["secondary_or_supplementary"] += 1
                continue
            if read.is_qcfail or read.is_duplicate:
                summary["qcfail_or_duplicate"] += 1
                continue
            if read.mapping_quality < args.min_mapq:
                summary["low_mapq"] += 1
                continue
            if read.query_sequence is None:
                summary["no_sequence"] += 1
                continue

            cb = get_cb(read)
            if cb is None:
                summary["missing_CB"] += 1
                continue
            if whitelist is not None and cb not in whitelist:
                summary["not_in_whitelist"] += 1
                continue

            reads_used += 1
            cb_reads[cb] += 1

            if args.max_reads and reads_used > args.max_reads:
                break

            pairs = read.get_aligned_pairs(matches_only=True)
            if not pairs:
                continue

            ref_positions = [rpos for qpos, rpos in pairs if rpos is not None]
            if not ref_positions:
                continue

            chrom = bam.get_reference_name(read.reference_id)
            if chrom not in bin_map:
                summary["chrom_not_in_bins"] += 1
                continue

            start = max(0, min(ref_positions) - 1)
            end = max(ref_positions) + 2

            try:
                refseq = fa.fetch(chrom, start, end).upper()
            except Exception:
                summary["fasta_fetch_fail"] += 1
                continue

            qseq = read.query_sequence.upper()
            qqual = read.query_qualities

            for qpos, rpos in pairs:
                if qpos is None or rpos is None:
                    continue

                if qqual is not None and qqual[qpos] < args.min_baseq:
                    continue

                off = rpos - start
                if off < 0 or off >= len(refseq):
                    continue

                ref_base = refseq[off]
                read_base = qseq[qpos]
                bin_idx = rpos // 100000
                bin_id = bin_map[chrom].get(bin_idx)

                if bin_id is None:
                    continue

                modified = None

                # CpG on reference C position: C->T as TAPS modified-like signal
                if ref_base == "C":
                    if off + 1 < len(refseq) and refseq[off + 1] == "G":
                        if read_base == "T":
                            modified = 1
                        elif read_base == "C":
                            modified = 0

                # CpG on opposite strand reference G position: G->A as modified-like signal
                elif ref_base == "G":
                    if off - 1 >= 0 and refseq[off - 1] == "C":
                        if read_base == "A":
                            modified = 1
                        elif read_base == "G":
                            modified = 0

                if modified is None:
                    continue

                cb_cov[cb] += 1
                cb_bin_cov[(cb, bin_id)] += 1

                if modified == 1:
                    cb_mod[cb] += 1
                    cb_bin_mod[(cb, bin_id)] += 1

            if args.max_reads and reads_used >= args.max_reads:
                break

        bam.close()

        if args.max_reads and reads_used >= args.max_reads:
            break

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    global_out = str(out_prefix) + ".global_per_cb.tsv"
    bin_out = str(out_prefix) + ".bin_counts.tsv.gz"
    summary_out = str(out_prefix) + ".summary.tsv"

    covered_bins_by_cb = Counter()
    for cb, bin_id in cb_bin_cov:
        covered_bins_by_cb[cb] += 1

    with open(global_out, "w") as out:
        out.write("CB\treads_used\tmodified_calls\tcovered_calls\tmodified_rate\tcovered_bins\n")
        for cb in sorted(cb_cov.keys(), key=lambda x: cb_cov[x], reverse=True):
            cov = cb_cov[cb]
            mod = cb_mod[cb]
            rate = mod / cov if cov else 0
            out.write(f"{cb}\t{cb_reads[cb]}\t{mod}\t{cov}\t{rate:.6f}\t{covered_bins_by_cb[cb]}\n")

    with gzip.open(bin_out, "wt", compresslevel=1) as out:
        out.write("CB\tbin_id\tmodified_calls\tcovered_calls\tmodified_rate\n")
        for key, cov in cb_bin_cov.items():
            cb, bin_id = key
            mod = cb_bin_mod[key]
            rate = mod / cov if cov else 0
            out.write(f"{cb}\t{bin_id}\t{mod}\t{cov}\t{rate:.6f}\n")

    with open(summary_out, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"bam_files\t{len(bam_paths)}\n")
        out.write(f"reads_used\t{reads_used}\n")
        out.write(f"cells_with_covered_CpG\t{len(cb_cov)}\n")
        total_mod = sum(cb_mod.values())
        total_cov = sum(cb_cov.values())
        out.write(f"total_modified_calls\t{total_mod}\n")
        out.write(f"total_covered_calls\t{total_cov}\n")
        out.write(f"overall_modified_rate\t{(total_mod / total_cov if total_cov else 0):.6f}\n")
        for k, v in summary.items():
            out.write(f"{k}\t{v}\n")

    sys.stderr.write(f"[INFO] wrote {global_out}\n")
    sys.stderr.write(f"[INFO] wrote {bin_out}\n")
    sys.stderr.write(f"[INFO] wrote {summary_out}\n")

if __name__ == "__main__":
    main()
