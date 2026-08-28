#!/usr/bin/env python3
import re
import pandas as pd
import pysam
from pathlib import Path

A_UMAP = "08_qc_summary/final_single_sample/region_umap_paper/P7_no01_promoter_residual_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"
BAM_LIST = "08_qc_summary/final_single_sample/A_cluster_bw_bam.list"

cb_set = set(pd.read_csv(A_UMAP, sep="\t")["CB"].astype(str))
print("[INFO] A CB count:", len(cb_set))
print("[INFO] A CB examples:", list(sorted(cb_set))[:5])

bams = [x.strip() for x in open(BAM_LIST) if x.strip() and not x.startswith("#")]
print("[INFO] bam count:", len(bams))

for bam_path in bams[:20]:
    print("\n==============================")
    print("[BAM]", bam_path)

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        print("[ERROR] cannot open:", e)
        continue

    print("[INFO] references head:", bam.references[:10])

    n = 0
    mapped = 0
    xm_tag = 0
    cb_tag_any = 0
    cb_tag_hit = 0
    prefix36_hit = 0
    any_36mer_hit = 0

    examples = []

    for read in bam.fetch(until_eof=True):
        n += 1
        if read.is_unmapped:
            continue
        mapped += 1

        qn = read.query_name
        tags = dict(read.get_tags())

        if "XM" in tags:
            xm_tag += 1

        tag_has_cb = False
        for tag in ["CB", "CR", "BX", "XC"]:
            if tag in tags:
                tag_has_cb = True
                if str(tags[tag]) in cb_set:
                    cb_tag_hit += 1
                    break
        if tag_has_cb:
            cb_tag_any += 1

        if len(qn) >= 36 and qn[:36] in cb_set:
            prefix36_hit += 1

        found_36mer = False
        for m in re.finditer(r"[ACGTN]{36}", qn):
            if m.group(0) in cb_set:
                found_36mer = True
                break
        if found_36mer:
            any_36mer_hit += 1

        if len(examples) < 5:
            examples.append({
                "qname": qn,
                "prefix36": qn[:36],
                "tags": read.get_tags()[:15],
                "has_XM": "XM" in tags,
                "CB_tags": {t: tags.get(t) for t in ["CB", "CR", "BX", "XC"] if t in tags},
            })

        if mapped >= 20000:
            break

    bam.close()

    print("[SUMMARY first mapped reads]")
    print("checked_reads:", n)
    print("mapped_checked:", mapped)
    print("XM_tag_count:", xm_tag)
    print("any_CB_like_tag_count:", cb_tag_any)
    print("CB_tag_hit_A:", cb_tag_hit)
    print("qname_prefix36_hit_A:", prefix36_hit)
    print("qname_any36mer_hit_A:", any_36mer_hit)

    print("[EXAMPLES]")
    for ex in examples:
        print("qname:", ex["qname"])
        print("prefix36:", ex["prefix36"])
        print("has_XM:", ex["has_XM"])
        print("CB_tags:", ex["CB_tags"])
        print("tags:", ex["tags"])
        print("---")
