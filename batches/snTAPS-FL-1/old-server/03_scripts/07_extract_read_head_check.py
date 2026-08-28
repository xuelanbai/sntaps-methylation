#!/usr/bin/env python3
import gzip
from pathlib import Path
import re

PROJECT = Path("/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation")
RAWDIR = PROJECT / "01_rawdata"
OUTDIR = PROJECT / "08_qc_summary/read_structure"
OUTDIR.mkdir(parents=True, exist_ok=True)

r1_files = sorted(RAWDIR.glob("*/*_1.fq.gz"))[:8]

print("[INFO] R1 files selected:")
for x in r1_files:
    print(x)
print(f"[INFO] R1 count: {len(r1_files)}")

if not r1_files:
    raise SystemExit("[ERROR] no R1 fastq found")

def read_n_pairs(r1, r2, n=10):
    with gzip.open(r1, "rt") as f1, gzip.open(r2, "rt") as f2:
        for _ in range(n):
            h1 = f1.readline().rstrip()
            s1 = f1.readline().rstrip()
            p1 = f1.readline().rstrip()
            q1 = f1.readline().rstrip()

            h2 = f2.readline().rstrip()
            s2 = f2.readline().rstrip()
            p2 = f2.readline().rstrip()
            q2 = f2.readline().rstrip()

            if not h1 or not h2:
                break

            yield h1, s1, h2, s2

for r1 in r1_files:
    r2 = Path(str(r1).replace("_1.fq.gz", "_2.fq.gz"))
    if not r2.exists():
        raise SystemExit(f"[ERROR] missing R2: {r2}")

    sample = r1.parent.name
    m = re.search(r"_L(\d+)_1\.fq\.gz$", r1.name)
    lane = f"L{m.group(1)}" if m else "unknownLane"

    out = OUTDIR / f"{sample}_{lane}.head_check.txt"

    with open(out, "w") as fo:
        fo.write(f"### {sample} {lane}\n")
        fo.write(f"R1={r1}\n")
        fo.write(f"R2={r2}\n\n")

        for i, (h1, s1, h2, s2) in enumerate(read_n_pairs(r1, r2, 10), start=1):
            idx = h1.split()[-1].split(":")[-1] if h1.split() else "NA"
            fo.write(f"PAIR_{i}\n")
            fo.write(f"HEADER_R1\t{h1}\n")
            fo.write(f"HEADER_R2\t{h2}\n")
            fo.write(f"HEADER_INDEX\t{idx}\n")
            fo.write(f"R1_first40\t{s1[:40]}\n")
            fo.write(f"R2_first40\t{s2[:40]}\n")
            fo.write(f"R1_first8\t{s1[:8]}\n")
            fo.write(f"R2_first8\t{s2[:8]}\n\n")

    print(f"[INFO] wrote {out}")
