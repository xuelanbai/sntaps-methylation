#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from collections import defaultdict
import subprocess
import shutil

import pandas as pd
import pysam


def extract_cb(read, cb_set, barcode_len=36):
    # 优先读tag
    for tag in ["CB", "CR", "BX", "XC"]:
        try:
            v = read.get_tag(tag)
            if v in cb_set:
                return v
        except KeyError:
            pass

    qn = read.query_name

    # 常见情况：read name开头就是36bp CB
    if len(qn) >= barcode_len:
        v = qn[:barcode_len]
        if v in cb_set:
            return v

    # 再尝试从read name里找36bp ACGTN片段
    for m in re.finditer(r"[ACGTN]{30,44}", qn):
        v = m.group(0)
        if v in cb_set:
            return v

    # 最后按分隔符拆
    for tok in re.split(r"[:|_\s/]+", qn):
        if tok in cb_set:
            return tok

    return None


def choose_bam_chrom(bam, chrom):
    refs = set(bam.references)
    if chrom in refs:
        return chrom
    if chrom.startswith("chr") and chrom[3:] in refs:
        return chrom[3:]
    if ("chr" + chrom) in refs:
        return "chr" + chrom
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umap", required=True)
    ap.add_argument("--bam-list", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mincov", type=int, default=5)
    ap.add_argument("--mapq", type=int, default=10)
    ap.add_argument("--barcode-len", type=int, default=36)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    umap = pd.read_csv(args.umap, sep="\t")
    if "CB" not in umap.columns or "cluster" not in umap.columns:
        raise ValueError(f"UMAP needs CB and cluster columns, got {umap.columns.tolist()}")

    umap["CB"] = umap["CB"].astype(str)
    umap["cluster"] = umap["cluster"].astype(int)

    cb_to_cluster = dict(zip(umap["CB"], umap["cluster"]))
    cb_set = set(cb_to_cluster)
    clusters = sorted(umap["cluster"].unique())

    bams = [x.strip() for x in open(args.bam_list) if x.strip() and not x.startswith("#")]

    print(f"[INFO] selected cells: {len(cb_set)}", flush=True)
    print(f"[INFO] clusters: {clusters}", flush=True)
    print(f"[INFO] target chrom: {args.chrom}", flush=True)
    print(f"[INFO] bam count: {len(bams)}", flush=True)

    # agg[cluster][position_0based] = [modified, covered]
    agg = {c: defaultdict(lambda: [0, 0]) for c in clusters}

    total_reads = 0
    used_reads = 0
    reads_with_xm = 0
    reads_with_cb = 0
    cpg_obs = 0

    for bam_path in bams:
        print(f"[INFO] open bam: {bam_path}", flush=True)
        bam = pysam.AlignmentFile(bam_path, "rb")
        bam_chrom = choose_bam_chrom(bam, args.chrom)

        if bam_chrom is None:
            print(f"[WARN] chrom {args.chrom} not found in {bam_path}, skip", flush=True)
            bam.close()
            continue

        for read in bam.fetch(bam_chrom):
            total_reads += 1

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < args.mapq:
                continue

            cb = extract_cb(read, cb_set, args.barcode_len)
            if cb is None:
                continue
            reads_with_cb += 1

            try:
                xm = read.get_tag("XM")
            except KeyError:
                continue
            reads_with_xm += 1

            cl = cb_to_cluster[cb]
            pairs = read.get_aligned_pairs(matches_only=False)

            for qpos, rpos in pairs:
                if qpos is None or rpos is None:
                    continue
                if qpos >= len(xm):
                    continue

                code = xm[qpos]
                # Bismark CpG context: Z methylated CpG, z unmethylated CpG
                if code == "Z":
                    agg[cl][rpos][0] += 1
                    agg[cl][rpos][1] += 1
                    cpg_obs += 1
                elif code == "z":
                    agg[cl][rpos][1] += 1
                    cpg_obs += 1

            used_reads += 1

            if total_reads % 1000000 == 0:
                print(
                    f"[INFO] {args.chrom} reads={total_reads} used={used_reads} cb={reads_with_cb} xm={reads_with_xm} cpg_obs={cpg_obs}",
                    flush=True
                )

        bam.close()

    print(f"[INFO] total_reads={total_reads}", flush=True)
    print(f"[INFO] reads_with_cb={reads_with_cb}", flush=True)
    print(f"[INFO] reads_with_xm={reads_with_xm}", flush=True)
    print(f"[INFO] used_reads={used_reads}", flush=True)
    print(f"[INFO] cpg_obs={cpg_obs}", flush=True)

    for cl in clusters:
        methyl_bg = outdir / f"C{cl}.{args.chrom}.methyl.bedGraph"
        hypo_bg = outdir / f"C{cl}.{args.chrom}.hypo.bedGraph"
        cov_bg = outdir / f"C{cl}.{args.chrom}.coverage.bedGraph"

        n = 0
        with open(methyl_bg, "w") as fm, open(hypo_bg, "w") as fh, open(cov_bg, "w") as fc:
            for pos in sorted(agg[cl]):
                mod, cov = agg[cl][pos]
                if cov < args.mincov:
                    continue
                meth = mod / cov
                hypo = 1.0 - meth

                # bedGraph: chrom start end value, 0-based half-open
                fm.write(f"{args.chrom}\t{pos}\t{pos+1}\t{meth:.6f}\n")
                fh.write(f"{args.chrom}\t{pos}\t{pos+1}\t{hypo:.6f}\n")
                fc.write(f"{args.chrom}\t{pos}\t{pos+1}\t{cov}\n")
                n += 1

        print(f"[INFO] C{cl} wrote {n} CpG sites on {args.chrom}", flush=True)


if __name__ == "__main__":
    main()
