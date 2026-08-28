#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import re

BASE = Path(".")
OUT = BASE / "08_qc_summary/final_single_sample/marker_pig_frontal_lobe/debug_marker_overlap"
OUT.mkdir(parents=True, exist_ok=True)

MARKER = BASE / "03_scripts/annotation/pig_frontal_lobe_marker_symbols.tsv"
BED = BASE / "02_reference/pig/gene_regions/pig_promoter_TSS_1.5kb.bed"

RAW = BASE / "08_qc_summary/final_single_sample/region_matrix/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.region_counts.tsv.gz"
COV10 = BASE / "08_qc_summary/final_single_sample/region_matrix_covfilter/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.cov10.region_counts.tsv.gz"

CELL = BASE / "08_qc_summary/final_single_sample/cell_filter_P7_QCpass/snTAPS-FL-1.P7_QCpass_strict.global_per_cb.tsv"

markers = pd.read_csv(MARKER, sep="\t")
markers["gene_upper"] = markers["gene"].astype(str).str.upper()
marker_set = set(markers["gene_upper"])

print("[INFO] marker genes:", len(marker_set))

# 解析BED，尽量同时记录region_id和gene symbol
bed_rows = []
with open(BED) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = re.split(r"\s+", line.strip())

        # 常见BED前4列：chrom start end name
        region_id = parts[3] if len(parts) >= 4 else None

        # 同时尝试找ENSSSCG
        ens = [x for x in parts if re.match(r"^ENSSSCG\d+", x)]
        if ens:
            region_id = ens[0]

        hits = []
        for x in parts:
            xu = x.upper()
            if xu in marker_set:
                hits.append(xu)

        for g in sorted(set(hits)):
            bed_rows.append((region_id, g, line.strip()))

bed_map = pd.DataFrame(bed_rows, columns=["region_id_from_bed", "gene_upper", "bed_line"]).drop_duplicates()
bed_map.to_csv(OUT / "marker_regions_from_bed.tsv", sep="\t", index=False)

print("[INFO] marker regions from BED:", bed_map.shape)
print("[INFO] found marker genes:", bed_map["gene_upper"].nunique() if not bed_map.empty else 0)

if bed_map.empty:
    print("[ERROR] no marker genes found in BED. BED may not contain gene symbols.")
    raise SystemExit(1)

marker_regions = set(bed_map["region_id_from_bed"].astype(str))

cells = pd.read_csv(CELL, sep="\t", usecols=["CB"])
cells["CB"] = cells["CB"].astype(str)
cell_set = set(cells["CB"])
print("[INFO] selected cells:", len(cell_set))

def inspect_matrix(path, label):
    print(f"\n[INFO] inspect {label}: {path}")
    head = pd.read_csv(path, sep="\t", nrows=5)
    print("[INFO] columns:", head.columns.tolist())
    print(head.head().to_string(index=False))

    region_col = "region_id" if "region_id" in head.columns else head.columns[1]

    all_regions = set()
    hit_regions = set()
    hit_rows_region_only = 0
    hit_rows_cell_only = 0
    hit_rows_both = 0

    n = 0
    for chunk in pd.read_csv(path, sep="\t", chunksize=1000000):
        n += len(chunk)
        chunk["CB"] = chunk["CB"].astype(str)
        chunk[region_col] = chunk[region_col].astype(str)

        regs = set(chunk[region_col].unique())
        all_regions.update(regs)
        hit_regions.update(regs & marker_regions)

        mr = chunk[region_col].isin(marker_regions)
        mc = chunk["CB"].isin(cell_set)

        hit_rows_region_only += int(mr.sum())
        hit_rows_cell_only += int(mc.sum())
        hit_rows_both += int((mr & mc).sum())

    print("[INFO]", label, "rows:", n)
    print("[INFO]", label, "unique regions:", len(all_regions))
    print("[INFO]", label, "marker regions overlapping matrix:", len(hit_regions))
    print("[INFO]", label, "rows with marker region:", hit_rows_region_only)
    print("[INFO]", label, "rows with selected cells:", hit_rows_cell_only)
    print("[INFO]", label, "rows with marker region AND selected cells:", hit_rows_both)

    pd.Series(sorted(list(hit_regions))).to_csv(OUT / f"{label}.overlap_marker_regions.txt", index=False, header=False)

inspect_matrix(RAW, "raw")
inspect_matrix(COV10, "cov10")
