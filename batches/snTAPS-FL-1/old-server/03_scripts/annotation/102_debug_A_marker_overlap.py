#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

A_UMAP = Path("08_qc_summary/final_single_sample/region_umap_paper/P7_no01_promoter_residual_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv")

REGION_MAP = Path("08_qc_summary/final_single_sample/marker_pig_frontal_lobe/A_no01_rerun_marker/A_no01_marker_region_map.tsv")

MATRIX_FULL = Path("08_qc_summary/final_single_sample/region_matrix/snTAPS-FL-1.strict2.promoter_TSS_1.5kb.region_counts.tsv.gz")
MATRIX_NO01 = Path("08_qc_summary/final_single_sample/region_matrix_P7_no01/snTAPS-FL-1.P7_no01.promoter_TSS_1.5kb.region_counts.tsv.gz")

QC = Path("08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.global_per_cb.tsv")

umap = pd.read_csv(A_UMAP, sep="\t")
regions = pd.read_csv(REGION_MAP, sep="\t")

print("### A_UMAP")
print("path:", A_UMAP)
print("columns:", umap.columns.tolist())
print("n_cells:", umap.shape[0])
print("cluster counts:")
print(umap["cluster"].value_counts().sort_index())
print("CB examples:")
print(umap["CB"].head(10).to_string(index=False))
print("CB length summary:")
print(umap["CB"].astype(str).str.len().describe())

print("\n### marker region map")
print("path:", REGION_MAP)
print("n_marker_regions:", regions["region_id"].nunique())
print("region examples:")
print(regions["region_id"].head(10).to_string(index=False))

# check CB overlap with global QC
qc = pd.read_csv(QC, sep="\t", usecols=["CB"])
a_cells = set(umap["CB"].astype(str))
qc_cells = set(qc["CB"].astype(str))
print("\n### CB overlap with global QC")
print("A cells:", len(a_cells))
print("QC cells:", len(qc_cells))
print("A ∩ QC:", len(a_cells & qc_cells))

marker_regions = set(regions["region_id"].astype(str))

def scan_matrix(path):
    print(f"\n### scan matrix: {path}")
    if not path.exists():
        print("NOT EXISTS")
        return

    total = 0
    cb_hit = 0
    region_hit = 0
    both_hit = 0
    example_both = []

    for i, chunk in enumerate(pd.read_csv(path, sep="\t", usecols=["CB", "region_id"], chunksize=1000000), 1):
        chunk["CB"] = chunk["CB"].astype(str)
        chunk["region_id"] = chunk["region_id"].astype(str)

        c = chunk["CB"].isin(a_cells)
        r = chunk["region_id"].isin(marker_regions)

        total += len(chunk)
        cb_hit += int(c.sum())
        region_hit += int(r.sum())
        both_hit += int((c & r).sum())

        if len(example_both) < 5:
            ex = chunk.loc[c & r].head(5)
            if not ex.empty:
                example_both.extend(ex.to_dict("records"))

        print(
            f"[chunk {i}] total={total} cb_hit={cb_hit} region_hit={region_hit} both_hit={both_hit}",
            flush=True
        )

    print("FINAL total_rows:", total)
    print("FINAL rows_with_A_cells:", cb_hit)
    print("FINAL rows_with_marker_regions:", region_hit)
    print("FINAL rows_with_A_cells_AND_marker_regions:", both_hit)
    print("example both rows:", example_both[:5])

scan_matrix(MATRIX_NO01)
scan_matrix(MATRIX_FULL)
