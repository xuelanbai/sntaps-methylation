#!/usr/bin/env python3
import argparse
from pathlib import Path
from collections import defaultdict

import pandas as pd
import pysam


def get_cb(read, cb_set):
    try:
        raw = str(read.get_tag("CB"))
        cb = raw.split()[0]
        if cb in cb_set:
            return cb
    except KeyError:
        return None
    return None


def choose_chrom(bam, chrom):
    refs = set(bam.references)
    if chrom in refs:
        return chrom
    if ("chr" + chrom) in refs:
        return "chr" + chrom
    if chrom.startswith("chr") and chrom[3:] in refs:
        return chrom[3:]
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umap", required=True)
    ap.add_argument("--bam-list", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mincov", type=int, default=5)
    ap.add_argument("--mapq", type=int, default=10)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    umap = pd.read_csv(args.umap, sep="\t")
    umap["CB"] = umap["CB"].astype(str)
    umap["cluster"] = pd.to_numeric(umap["cluster"], errors="coerce").astype(int)

    cb_to_cluster = dict(zip(umap["CB"], umap["cluster"]))
    cb_set = set(cb_to_cluster)
    clusters = sorted(umap["cluster"].unique())

    bams = [x.strip() for x in open(args.bam_list) if x.strip() and not x.startswith("#")]

    fa = pysam.FastaFile(args.fasta)
    if args.chrom not in fa.references and ("chr" + args.chrom) not in fa.references:
        raise ValueError(f"chrom {args.chrom} not found in fasta. First refs={fa.references[:10]}")

    fa_chrom = args.chrom if args.chrom in fa.references else "chr" + args.chrom
    ref = fa.fetch(fa_chrom).upper()
    chrom_len = len(ref)

    print(f"[INFO] selected cells: {len(cb_set)}", flush=True)
    print(f"[INFO] clusters: {clusters}", flush=True)
    print(f"[INFO] chrom: {args.chrom} fasta_chrom: {fa_chrom} length={chrom_len}", flush=True)
    print(f"[INFO] bam count: {len(bams)}", flush=True)

    # agg[cluster][CpG_C_0based] = [modified, covered]
    agg = {c: defaultdict(lambda: [0, 0]) for c in clusters}

    total_reads = 0
    mapped_reads = 0
    cb_reads = 0
    used_reads = 0
    cpg_obs = 0

    for bam_path in bams:
        print(f"[INFO] open bam: {bam_path}", flush=True)
        bam = pysam.AlignmentFile(bam_path, "rb")
        bam_chrom = choose_chrom(bam, args.chrom)
        if bam_chrom is None:
            print(f"[WARN] chrom {args.chrom} not in bam, skip {bam_path}", flush=True)
            bam.close()
            continue

        for read in bam.fetch(bam_chrom):
            total_reads += 1

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < args.mapq:
                continue

            mapped_reads += 1

            cb = get_cb(read, cb_set)
            if cb is None:
                continue
            cb_reads += 1

            qseq = read.query_sequence
            if not qseq:
                continue
            qseq = qseq.upper()

            cl = cb_to_cluster[cb]
            read_used = False

            for qpos, rpos in read.get_aligned_pairs(matches_only=False):
                if qpos is None or rpos is None:
                    continue
                if rpos < 0 or rpos >= chrom_len:
                    continue
                if qpos < 0 or qpos >= len(qseq):
                    continue

                qb = qseq[qpos]
                rb = ref[rpos]

                # Plus-strand CpG C: reference C followed by G.
                # C = unmodified/covered, T = modified/covered.
                if rb == "C" and rpos + 1 < chrom_len and ref[rpos + 1] == "G":
                    if qb == "C":
                        agg[cl][rpos][1] += 1
                        cpg_obs += 1
                        read_used = True
                    elif qb == "T":
                        agg[cl][rpos][0] += 1
                        agg[cl][rpos][1] += 1
                        cpg_obs += 1
                        read_used = True

                # Reverse-complement CpG C corresponds to reference G preceded by C.
                # G = unmodified/covered, A = modified/covered.
                elif rb == "G" and rpos - 1 >= 0 and ref[rpos - 1] == "C":
                    cpg_c_pos = rpos - 1
                    if qb == "G":
                        agg[cl][cpg_c_pos][1] += 1
                        cpg_obs += 1
                        read_used = True
                    elif qb == "A":
                        agg[cl][cpg_c_pos][0] += 1
                        agg[cl][cpg_c_pos][1] += 1
                        cpg_obs += 1
                        read_used = True

            if read_used:
                used_reads += 1

            if total_reads % 1000000 == 0:
                print(
                    f"[INFO] reads={total_reads} mapped={mapped_reads} cb={cb_reads} used={used_reads} cpg_obs={cpg_obs}",
                    flush=True
                )

        bam.close()

    print(f"[INFO] total_reads={total_reads}", flush=True)
    print(f"[INFO] mapped_reads={mapped_reads}", flush=True)
    print(f"[INFO] cb_reads={cb_reads}", flush=True)
    print(f"[INFO] used_reads={used_reads}", flush=True)
    print(f"[INFO] cpg_obs={cpg_obs}", flush=True)

    for cl in clusters:
        methyl_bg = outdir / f"C{cl}.{args.chrom}.methyl.bedGraph"
        hypo_bg = outdir / f"C{cl}.{args.chrom}.hypo.bedGraph"
        cov_bg = outdir / f"C{cl}.{args.chrom}.coverage.bedGraph"

        n = 0
        with open(methyl_bg, "w") as fm, open(hypo_bg, "w") as fh, open(cov_bg, "w") as fc:
            for pos in sorted(agg[cl].keys()):
                mod, cov = agg[cl][pos]
                if cov < args.mincov:
                    continue
                meth = mod / cov
                hypo = 1.0 - meth

                fm.write(f"{args.chrom}\t{pos}\t{pos+1}\t{meth:.6f}\n")
                fh.write(f"{args.chrom}\t{pos}\t{pos+1}\t{hypo:.6f}\n")
                fc.write(f"{args.chrom}\t{pos}\t{pos+1}\t{cov}\n")
                n += 1

        print(f"[INFO] C{cl} wrote {n} CpG sites on {args.chrom}", flush=True)

    fa.close()


if __name__ == "__main__":
    main()
