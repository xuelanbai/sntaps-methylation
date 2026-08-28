#!/usr/bin/env python3
import gzip
import sys
from collections import Counter

def load_barcode(path):
    d = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            p = line.split()
            if len(p) >= 2:
                name = p[0]
                seq = p[1].upper()
                d[seq] = name
    return d

def fq_iter(r1, r2, n):
    with gzip.open(r1, "rt") as f1, gzip.open(r2, "rt") as f2:
        for _ in range(n):
            h1 = f1.readline().strip()
            if not h1:
                break
            s1 = f1.readline().strip().upper()
            f1.readline()
            f1.readline()

            h2 = f2.readline().strip()
            s2 = f2.readline().strip().upper()
            f2.readline()
            f2.readline()

            yield h1, s1, h2, s2

def get_header_index(header):
    # 取最后一个冒号后的 index，例如 CACTAGCCTT+TGGCAGACTA
    return header.split()[-1].split(":")[-1].upper()

def main():
    if len(sys.argv) < 7:
        print("Usage: python3 script.py R1.fq.gz R2.fq.gz T5 T7 I5 I7 [N]", file=sys.stderr)
        sys.exit(1)

    r1, r2, t5_file, t7_file, i5_file, i7_file = sys.argv[1:7]
    n = int(sys.argv[7]) if len(sys.argv) > 7 else 100000

    T5 = load_barcode(t5_file)
    T7 = load_barcode(t7_file)
    I5 = load_barcode(i5_file)
    I7 = load_barcode(i7_file)

    c = Counter()
    top_r1 = Counter()
    top_r2 = Counter()
    top_idx = Counter()

    for h1, s1, h2, s2 in fq_iter(r1, r2, n):
        c["total"] += 1

        r1_8 = s1[:8]
        r2_8 = s2[:8]
        idx = get_header_index(h1)
        parts = idx.split("+")

        idx1 = parts[0] if len(parts) > 0 else ""
        idx2 = parts[1] if len(parts) > 1 else ""

        top_r1[r1_8] += 1
        top_r2[r2_8] += 1
        top_idx[idx] += 1

        if r1_8 in T5:
            c["R1_first8_in_T5"] += 1
        if r1_8 in T7:
            c["R1_first8_in_T7"] += 1
        if r2_8 in T5:
            c["R2_first8_in_T5"] += 1
        if r2_8 in T7:
            c["R2_first8_in_T7"] += 1

        if idx1 in I7 and idx2 in I5:
            c["idx1_I7_idx2_I5"] += 1
        if idx1 in I5 and idx2 in I7:
            c["idx1_I5_idx2_I7"] += 1

        # 旧scATAC假设：R2前8=T7，R1前8=T5，header=index1(I7)+index2(I5)
        if r2_8 in T7 and r1_8 in T5 and idx1 in I7 and idx2 in I5:
            c["old_scATAC_layout_match"] += 1

    total = c["total"] if c["total"] else 1

    print(f"# R1\t{r1}")
    print(f"# R2\t{r2}")
    print(f"# checked_pairs\t{c['total']}")
    print("metric\tcount\tpercent")

    for k in [
        "R1_first8_in_T5",
        "R1_first8_in_T7",
        "R2_first8_in_T5",
        "R2_first8_in_T7",
        "idx1_I7_idx2_I5",
        "idx1_I5_idx2_I7",
        "old_scATAC_layout_match",
    ]:
        print(f"{k}\t{c[k]}\t{c[k] / total * 100:.2f}")

    print("\n# top_R1_first8")
    for seq, cnt in top_r1.most_common(20):
        print(seq, cnt, T5.get(seq, "-"), T7.get(seq, "-"), sep="\t")

    print("\n# top_R2_first8")
    for seq, cnt in top_r2.most_common(20):
        print(seq, cnt, T5.get(seq, "-"), T7.get(seq, "-"), sep="\t")

    print("\n# top_header_index")
    for seq, cnt in top_idx.most_common(20):
        print(seq, cnt, sep="\t")

if __name__ == "__main__":
    main()
