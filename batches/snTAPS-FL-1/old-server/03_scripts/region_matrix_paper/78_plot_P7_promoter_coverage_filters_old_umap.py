#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(".")
P7 = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_residual_all_dropPC1_smooth100/snTAPS-FL-1.P7_promoter_residual_all_dropPC1_smooth100.umap.tsv"
CELL_DIR = BASE / "08_qc_summary/final_single_sample/cell_filter_P7_promoter_coverage"
OUT = BASE / "08_qc_summary/final_single_sample/region_umap_paper/P7_promoter_coverage_filters_old_umap"
OUT.mkdir(parents=True, exist_ok=True)

u = pd.read_csv(P7, sep="\t")
u["CB"] = u["CB"].astype(str)

for f in sorted(CELL_DIR.glob("snTAPS-FL-1.P7_*.global_per_cb.tsv")):
    name = f.name.replace("snTAPS-FL-1.", "").replace(".global_per_cb.tsv", "")
    cells = pd.read_csv(f, sep="\t")
    cells["CB"] = cells["CB"].astype(str)

    d = u[u["CB"].isin(set(cells["CB"]))].copy()
    print("[INFO]", name, len(d), flush=True)

    plt.figure(figsize=(7,6))
    sc = plt.scatter(d["UMAP1"], d["UMAP2"], c=d["cluster"], s=3, cmap="tab10")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(f"{name} on old P7 UMAP")
    plt.colorbar(sc, label="old P7 cluster")
    plt.tight_layout()
    plt.savefig(OUT / f"{name}.oldP7_by_cluster.png", dpi=220)
    plt.close()

    d["cluster"].value_counts().sort_index().to_csv(
        OUT / f"{name}.oldP7_cluster_counts.tsv",
        sep="\t",
        header=["n_cells"]
    )

print("[INFO] done", flush=True)
