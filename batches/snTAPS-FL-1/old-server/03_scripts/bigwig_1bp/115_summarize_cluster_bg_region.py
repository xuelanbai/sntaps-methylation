#!/usr/bin/env python3
import sys
from pathlib import Path

bg_dir = Path(sys.argv[1])
chrom = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
region = sys.argv[5]
out = Path(sys.argv[6])

rows = []
for c in range(6):
    meth_file = bg_dir / f"C{c}.methyl.sorted.bedGraph"
    cov_file = bg_dir / f"C{c}.coverage.sorted.bedGraph"

    meth = {}
    with open(meth_file) as f:
        for line in f:
            a = line.rstrip("\n").split("\t")
            if len(a) < 4:
                continue
            if a[0] != chrom:
                continue
            s = int(a[1]); e = int(a[2])
            if e <= start or s >= end:
                continue
            meth[(a[0], s, e)] = float(a[3])

    covered_sites = 0
    total_cov = 0.0
    weighted_mod = 0.0

    with open(cov_file) as f:
        for line in f:
            a = line.rstrip("\n").split("\t")
            if len(a) < 4:
                continue
            if a[0] != chrom:
                continue
            s = int(a[1]); e = int(a[2])
            if e <= start or s >= end:
                continue

            key = (a[0], s, e)
            if key not in meth:
                continue

            cov = float(a[3])
            m = meth[key]

            covered_sites += 1
            total_cov += cov
            weighted_mod += m * cov

    if total_cov > 0:
        weighted_methyl = weighted_mod / total_cov
        weighted_hypo = 1 - weighted_methyl
    else:
        weighted_methyl = "NA"
        weighted_hypo = "NA"

    rows.append([
        f"C{c}", region, chrom, start, end,
        covered_sites, total_cov, weighted_methyl, weighted_hypo
    ])

with open(out, "w") as w:
    w.write("cluster\tregion\tchrom\tstart\tend\tcovered_sites\ttotal_coverage\tweighted_methyl\tweighted_hypo\n")
    for r in rows:
        w.write("\t".join(map(str, r)) + "\n")
