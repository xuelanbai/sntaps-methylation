#!/usr/bin/env python3
import os
import glob
import gzip
import json
import subprocess
from pathlib import Path
from collections import defaultdict

PROJECT = Path("/share/org/YZWL/yzwl_wujx/Projects/snTAPS_methylation")
OUTDIR = PROJECT / "08_qc_summary/stage_summary"
OUTDIR.mkdir(parents=True, exist_ok=True)

samples = [f"snTAPS-FL-1_snTAPS-FL-1-{i}" for i in range(1, 9)]

def read_summary_tsv(path):
    d = {}
    if not path.exists():
        return d
    with open(path) as f:
        header = next(f, None)
        for line in f:
            arr = line.rstrip("\n").split("\t")
            if len(arr) >= 2:
                d[arr[0]] = arr[1]
    return d

def count_lines(path):
    if not path.exists():
        return "NA"
    n = 0
    with open(path) as f:
        for _ in f:
            n += 1
    return n

def safe_float(x):
    try:
        return float(x)
    except Exception:
        return None

def safe_int(x):
    try:
        return int(float(x))
    except Exception:
        return None

# 1) split stats
# expected rows like:
# snTAPS-FL-1_snTAPS-FL-1-1_L1.split.stat.tsv 281310840 85.78
split_file = PROJECT / "08_qc_summary/split_assigned_rate.tsv"
split_stats = defaultdict(lambda: {"lanes":0, "assigned":0, "total_est":0.0})

if split_file.exists():
    with open(split_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("sample"):
                continue
            arr = line.split()
            if len(arr) < 3:
                continue
            name, assigned, pct = arr[0], arr[1], arr[2]
            # remove .split.stat.tsv
            name = name.replace(".split.stat.tsv", "")
            # sample is before _L
            if "_L" not in name:
                continue
            sample = name.rsplit("_L", 1)[0]
            try:
                assigned = int(float(assigned))
                pct = float(pct)
                total_est = assigned / (pct / 100.0) if pct > 0 else 0
            except Exception:
                continue
            split_stats[sample]["lanes"] += 1
            split_stats[sample]["assigned"] += assigned
            split_stats[sample]["total_est"] += total_est

# 2) mapping stats from bam idxstats
map_stats = {}
for sample in samples:
    bam_dir = PROJECT / f"06_align/bam_by_lane/{sample}"
    bams = sorted(glob.glob(str(bam_dir / "*.sorted.bam")))
    mapped_total = 0
    unmapped_total = 0
    ok_lanes = 0

    for bam in bams:
        try:
            p = subprocess.run(
                ["samtools", "idxstats", bam],
                text=True,
                capture_output=True,
                check=True
            )
        except Exception:
            continue

        ok_lanes += 1
        for line in p.stdout.splitlines():
            arr = line.split("\t")
            if len(arr) < 4:
                continue
            mapped_total += int(arr[2])
            unmapped_total += int(arr[3])

    denom = mapped_total + unmapped_total
    mapped_pct = mapped_total / denom * 100 if denom > 0 else None
    map_stats[sample] = {
        "bam_lanes": ok_lanes,
        "mapped_reads": mapped_total,
        "unmapped_reads": unmapped_total,
        "mapped_pct": mapped_pct,
    }

# 3) methylation / merged summary and cell-level QC
rows = []

for sample in samples:
    outdir = PROJECT / f"08_qc_summary/downstream_test/{sample}"
    merged = outdir / "merged_full"

    wl = outdir / f"{sample}.CB.whitelist.readpairs_ge1000.txt"
    whitelist_cb = count_lines(wl)

    summary = merged / f"{sample}.full_merged.summary.tsv"
    global_cb = merged / f"{sample}.full_merged.global_per_cb.tsv"
    bin_counts = merged / f"{sample}.full_merged.bin_counts.tsv.gz"

    summ = read_summary_tsv(summary)

    # high-quality cells for clustering
    hq_cells = "NA"
    hq_pct = "NA"
    if global_cb.exists():
        total_cells = 0
        hq = 0
        with open(global_cb) as f:
            header = next(f, None)
            for line in f:
                arr = line.rstrip("\n").split("\t")
                if len(arr) < 6:
                    continue
                total_cells += 1
                # columns:
                # CB reads_used modified_calls covered_calls modified_rate covered_bins
                covered_calls = safe_int(arr[3])
                covered_bins = safe_int(arr[5])
                if covered_calls is not None and covered_bins is not None:
                    if covered_calls >= 10000 and covered_bins >= 500:
                        hq += 1
        hq_cells = hq
        hq_pct = hq / total_cells * 100 if total_cells > 0 else "NA"

    # merged bin count lines, optional but useful to detect incomplete merge
    bin_rows = "NA"
    if bin_counts.exists():
        try:
            p = subprocess.run(
                f"zcat {bin_counts} | wc -l",
                shell=True,
                text=True,
                capture_output=True,
                check=True
            )
            bin_rows = int(p.stdout.strip())
        except Exception:
            bin_rows = "NA"

    split = split_stats[sample]
    split_pct = (
        split["assigned"] / split["total_est"] * 100
        if split["total_est"] > 0 else None
    )

    m = map_stats.get(sample, {})

    rows.append({
        "sample": sample,
        "split_lanes": split["lanes"] if split["lanes"] else "NA",
        "split_total_pairs_est": int(split["total_est"]) if split["total_est"] else "NA",
        "split_assigned_pairs": split["assigned"] if split["assigned"] else "NA",
        "split_assigned_pct": f"{split_pct:.4f}" if split_pct is not None else "NA",
        "bam_lanes": m.get("bam_lanes", "NA"),
        "mapped_reads": m.get("mapped_reads", "NA"),
        "unmapped_reads": m.get("unmapped_reads", "NA"),
        "mapped_pct": f"{m['mapped_pct']:.4f}" if m.get("mapped_pct") is not None else "NA",
        "whitelist_cb_ge1000": whitelist_cb,
        "cells_with_covered_CpG": summ.get("cells_with_covered_CpG", "NA"),
        "hq_cells_cov10k_bins500": hq_cells,
        "hq_cell_pct": f"{hq_pct:.4f}" if isinstance(hq_pct, float) else "NA",
        "reads_used_in_calling": summ.get("reads_used", "NA"),
        "total_modified_calls": summ.get("total_modified_calls", "NA"),
        "total_covered_calls": summ.get("total_covered_calls", "NA"),
        "overall_modified_rate": summ.get("overall_modified_rate", "NA"),
        "merged_bin_count_rows": bin_rows,
    })

out = OUTDIR / "snTAPS_stage_summary_by_sample.tsv"
cols = [
    "sample",
    "split_lanes",
    "split_total_pairs_est",
    "split_assigned_pairs",
    "split_assigned_pct",
    "bam_lanes",
    "mapped_reads",
    "unmapped_reads",
    "mapped_pct",
    "whitelist_cb_ge1000",
    "cells_with_covered_CpG",
    "hq_cells_cov10k_bins500",
    "hq_cell_pct",
    "reads_used_in_calling",
    "total_modified_calls",
    "total_covered_calls",
    "overall_modified_rate",
    "merged_bin_count_rows",
]

with open(out, "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")

print(f"[INFO] wrote {out}")
