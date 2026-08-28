#!/usr/bin/env python3
from pathlib import Path
import argparse, re, shutil, subprocess
import pandas as pd
import numpy as np

def safe_name(x):
    x = re.sub(r"[^A-Za-z0-9]+", "_", x).strip("_")
    return x

def read_chrom_sizes(path):
    d = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            a = line.rstrip("\n").split("\t")
            d[a[0]] = int(a[1])
    return d

def collect_region_bed(region_bed, target_region_ids):
    records = {}
    with open(region_bed) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            a = re.split(r"\s+", line.strip())
            if len(a) < 3:
                continue
            chrom, start, end = a[0], int(a[1]), int(a[2])
            text = " ".join(a)
            ids = re.findall(r"ENSSSCG[0-9]+", text)

            # 优先选出现在matrix里的region_id
            hit = [x for x in ids if x in target_region_ids]
            if not hit:
                # 如果第4列本身就是region_id
                if len(a) >= 4 and a[3] in target_region_ids:
                    hit = [a[3]]
            for rid in hit:
                records[rid] = (chrom, start, end)
    return records

def load_marker_map(path):
    if path is None or not Path(path).exists():
        return None
    m = pd.read_csv(path, sep="\t")
    if "region_id" not in m.columns:
        return None
    if "label" not in m.columns:
        if {"group", "gene"}.issubset(m.columns):
            m["label"] = m["group"].astype(str) + ":" + m["gene"].astype(str)
        else:
            m["label"] = m["region_id"].astype(str)
    return m[["region_id", "label"]].drop_duplicates()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umap", required=True)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--region-bed", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--marker-map", default="")
    ap.add_argument("--min-cov", type=int, default=20)
    ap.add_argument("--make-bw", action="store_true")
    ap.add_argument("--bedGraphToBigWig", default="")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cluster_to_celltype = {
        0: "Inhibitory neurons",
        1: "Excitatory neurons I",
        2: "Vascular-associated cells",
        3: "Oligodendrocyte precursor cells",
        4: "Astrocytes",
        5: "Excitatory neurons II",
    }

    print("[INFO] read A UMAP", flush=True)
    umap = pd.read_csv(args.umap, sep="\t")
    umap["CB"] = umap["CB"].astype(str)
    umap["cluster"] = pd.to_numeric(umap["cluster"], errors="coerce").astype(int)

    cell2cluster = dict(zip(umap["CB"], umap["cluster"]))
    cells = set(cell2cluster)

    counts = umap["cluster"].value_counts().sort_index().rename_axis("cluster").reset_index(name="n_cells")
    counts["cell_type"] = counts["cluster"].map(cluster_to_celltype)
    counts.to_csv(outdir / "A_cluster_cell_counts.tsv", sep="\t", index=False)
    print(counts.to_string(index=False), flush=True)

    print("[INFO] scan matrix and aggregate promoter counts by cluster", flush=True)
    agg = {}
    usecols = ["CB", "region_id", "modified_calls", "covered_calls"]

    total = 0
    kept = 0

    for i, chunk in enumerate(pd.read_csv(args.matrix, sep="\t", usecols=usecols, chunksize=1000000), 1):
        chunk["CB"] = chunk["CB"].astype(str)
        chunk = chunk[chunk["CB"].isin(cells)].copy()
        if chunk.empty:
            continue

        chunk["cluster"] = chunk["CB"].map(cell2cluster)
        chunk["modified_calls"] = pd.to_numeric(chunk["modified_calls"], errors="coerce").fillna(0)
        chunk["covered_calls"] = pd.to_numeric(chunk["covered_calls"], errors="coerce").fillna(0)

        g = chunk.groupby(["cluster", "region_id"], as_index=False)[["modified_calls", "covered_calls"]].sum()

        for _, r in g.iterrows():
            key = (int(r["cluster"]), str(r["region_id"]))
            if key not in agg:
                agg[key] = [0.0, 0.0]
            agg[key][0] += float(r["modified_calls"])
            agg[key][1] += float(r["covered_calls"])

        total += len(chunk)
        kept += len(g)

        if i % 10 == 0:
            print(f"[INFO] chunks={i}, rows_from_A_cells={total}, aggregated_pairs_seen={kept}", flush=True)

    rows = []
    for (cl, rid), (mod, cov) in agg.items():
        if cov < args.min_cov:
            continue
        rate = mod / cov if cov > 0 else np.nan
        if np.isnan(rate):
            continue
        rows.append({
            "cluster": cl,
            "cell_type": cluster_to_celltype.get(cl, f"Cluster{cl}"),
            "region_id": rid,
            "modified_calls": mod,
            "covered_calls": cov,
            "methylation_rate": rate,
            "promoter_activity": 1.0 - rate,
        })

    merged = pd.DataFrame(rows)
    merged.to_csv(outdir / "A_cluster_promoter_merged_counts.tsv.gz", sep="\t", index=False, compression="gzip")
    print("[INFO] merged rows after min_cov:", merged.shape[0], flush=True)

    if merged.empty:
        raise RuntimeError("No merged promoter records. Check matrix / CB overlap / min_cov.")

    print("[INFO] map region_id to BED coordinates", flush=True)
    target_ids = set(merged["region_id"].astype(str))
    bedmap = collect_region_bed(args.region_bed, target_ids)
    print("[INFO] region IDs in merged:", len(target_ids), flush=True)
    print("[INFO] region IDs mapped to BED:", len(bedmap), flush=True)

    if len(bedmap) == 0:
        raise RuntimeError("No region_id mapped to BED. Need inspect promoter BED region_id column.")

    chrom_sizes = read_chrom_sizes(args.chrom_sizes)
    bw_tool = args.bedGraphToBigWig or shutil.which("bedGraphToBigWig")
    if args.make_bw and not bw_tool:
        print("[WARN] bedGraphToBigWig not found. Only bedGraph will be written.", flush=True)

    # marker promoter BED for IGV
    marker_map = load_marker_map(args.marker_map)
    if marker_map is not None:
        marker_rows = []
        for _, r in marker_map.iterrows():
            rid = str(r["region_id"])
            if rid in bedmap:
                chrom, start, end = bedmap[rid]
                marker_rows.append((chrom, start, end, str(r["label"])))
        marker_bed = outdir / "marker_promoters_for_IGV.bed"
        pd.DataFrame(marker_rows).drop_duplicates().to_csv(
            marker_bed, sep="\t", header=False, index=False
        )
        print("[INFO] wrote marker BED:", marker_bed, "n=", len(marker_rows), flush=True)

    print("[INFO] write bedGraph / bigWig per cluster", flush=True)

    for cl in sorted(merged["cluster"].unique()):
        celltype = cluster_to_celltype.get(int(cl), f"Cluster{cl}")
        name = f"C{cl}_{safe_name(celltype)}"
        sub = merged[merged["cluster"] == cl].copy()

        bed_rows = []
        missing = 0
        for _, r in sub.iterrows():
            rid = str(r["region_id"])
            if rid not in bedmap:
                missing += 1
                continue
            chrom, start, end = bedmap[rid]
            if chrom not in chrom_sizes:
                continue
            start = max(0, int(start))
            end = min(int(end), chrom_sizes[chrom])
            if end <= start:
                continue
            bed_rows.append((chrom, start, end, r["methylation_rate"], r["promoter_activity"]))

        bed = pd.DataFrame(bed_rows, columns=["chrom", "start", "end", "methylation_rate", "promoter_activity"])
        if bed.empty:
            print(f"[WARN] cluster {cl} has no BED rows after mapping", flush=True)
            continue

        # sort by chrom order in chrom.sizes then start
        chrom_order = {c:i for i,c in enumerate(chrom_sizes.keys())}
        bed["_ord"] = bed["chrom"].map(chrom_order)
        bed = bed.dropna(subset=["_ord"]).copy()
        bed["_ord"] = bed["_ord"].astype(int)
        bed = bed.sort_values(["_ord", "start", "end"])

        for track, col in [("methylation_rate", "methylation_rate"), ("promoter_activity", "promoter_activity")]:
            bg = outdir / f"{name}.{track}.bedGraph"
            bw = outdir / f"{name}.{track}.bw"

            bed[["chrom", "start", "end", col]].to_csv(bg, sep="\t", header=False, index=False)

            if args.make_bw and bw_tool:
                cmd = [bw_tool, str(bg), str(args.chrom_sizes), str(bw)]
                print("[INFO]", " ".join(cmd), flush=True)
                subprocess.run(cmd, check=True)

        print(f"[INFO] cluster {cl} {celltype}: BED rows={bed.shape[0]}, missing_region_bed={missing}", flush=True)

    print("[INFO] done. outdir:", outdir, flush=True)

if __name__ == "__main__":
    main()
