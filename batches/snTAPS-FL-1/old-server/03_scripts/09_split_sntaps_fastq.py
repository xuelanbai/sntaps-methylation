#!/usr/bin/env python3
import argparse
import gzip
from collections import Counter
from pathlib import Path

def load_bc(path):
    seq_to_name = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            p = line.split()
            if len(p) >= 2:
                name, seq = p[0], p[1].upper()
                seq_to_name[seq] = name
    return seq_to_name

def fq_iter(r1, r2):
    with gzip.open(r1, "rt") as f1, gzip.open(r2, "rt") as f2:
        while True:
            h1 = f1.readline().rstrip()
            if not h1:
                break
            s1 = f1.readline().rstrip().upper()
            p1 = f1.readline().rstrip()
            q1 = f1.readline().rstrip()

            h2 = f2.readline().rstrip()
            s2 = f2.readline().rstrip().upper()
            p2 = f2.readline().rstrip()
            q2 = f2.readline().rstrip()

            yield h1, s1, p1, q1, h2, s2, p2, q2

def header_index(header):
    return header.split()[-1].split(":")[-1].upper()

def write_read(out, header, seq, plus, qual):
    out.write(header + "\n")
    out.write(seq + "\n")
    out.write(plus + "\n")
    out.write(qual + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--out-r1", required=True)
    ap.add_argument("--out-r2", required=True)
    ap.add_argument("--t5", required=True)
    ap.add_argument("--t7", required=True)
    ap.add_argument("--i5", required=True)
    ap.add_argument("--i7", required=True)
    ap.add_argument("--trim-left", type=int, default=30)
    ap.add_argument("--max-reads", type=int, default=0)
    ap.add_argument("--stat", required=True)
    ap.add_argument("--bc-count", required=True)
    args = ap.parse_args()

    T5 = load_bc(args.t5)
    T7 = load_bc(args.t7)
    I5 = load_bc(args.i5)
    I7 = load_bc(args.i7)

    stats = Counter()
    bc_counter = Counter()

    Path(args.out_r1).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_r2).parent.mkdir(parents=True, exist_ok=True)
    Path(args.stat).parent.mkdir(parents=True, exist_ok=True)
    Path(args.bc_count).parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(args.out_r1, "wt", compresslevel=1) as o1, gzip.open(args.out_r2, "wt", compresslevel=1) as o2:
        for h1, s1, p1, q1, h2, s2, p2, q2 in fq_iter(args.r1, args.r2):
            stats["total_pairs"] += 1

            if args.max_reads and stats["total_pairs"] > args.max_reads:
                stats["total_pairs"] -= 1
                break

            if len(s1) <= args.trim_left or len(s2) <= args.trim_left:
                stats["too_short"] += 1
                continue

            r1_8 = s1[:8]
            r2_8 = s2[:8]
            idx = header_index(h1)
            parts = idx.split("+")

            if len(parts) != 2:
                stats["bad_header_index"] += 1
                continue

            idx1, idx2 = parts[0], parts[1]

            t5_name = T5.get(r1_8)
            t7_name = T7.get(r2_8)
            i7_name = I7.get(idx1)
            i5_name = I5.get(idx2)

            if t5_name is None:
                stats["fail_T5"] += 1
            if t7_name is None:
                stats["fail_T7"] += 1
            if i7_name is None:
                stats["fail_I7"] += 1
            if i5_name is None:
                stats["fail_I5"] += 1

            if not (t5_name and t7_name and i7_name and i5_name):
                stats["unassigned"] += 1
                continue

            # 旧scATAC逻辑：full barcode = T7 + T5 + I7 + I5
            cb_seq = r2_8 + r1_8 + idx1 + idx2
            cb_name = f"{t7_name}_{t5_name}_{i7_name}_{i5_name}"
            bc_counter[cb_name + "\t" + cb_seq] += 1

            new_h1 = h1 + f" CB:Z:{cb_seq} CN:Z:{cb_name}"
            new_h2 = h2 + f" CB:Z:{cb_seq} CN:Z:{cb_name}"

            new_s1 = s1[args.trim_left:]
            new_q1 = q1[args.trim_left:]
            new_s2 = s2[args.trim_left:]
            new_q2 = q2[args.trim_left:]

            write_read(o1, new_h1, new_s1, p1, new_q1)
            write_read(o2, new_h2, new_s2, p2, new_q2)

            stats["assigned"] += 1

    with open(args.stat, "w") as f:
        f.write("metric\tcount\tpercent\n")
        total = stats["total_pairs"] or 1
        for k in ["total_pairs", "assigned", "unassigned", "too_short", "bad_header_index", "fail_T5", "fail_T7", "fail_I7", "fail_I5"]:
            f.write(f"{k}\t{stats[k]}\t{stats[k]/total*100:.2f}\n")

    with open(args.bc_count, "w") as f:
        f.write("barcode_name\tbarcode_seq\tcount\n")
        for key, count in bc_counter.most_common():
            name, seq = key.split("\t")
            f.write(f"{name}\t{seq}\t{count}\n")

if __name__ == "__main__":
    main()
