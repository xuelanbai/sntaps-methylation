#!/usr/bin/env python3
import argparse
import gzip
import re
from pathlib import Path

def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def load_chrom_sizes(path):
    sizes = {}
    order = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            chrom, size = line.strip().split()[:2]
            sizes[chrom] = int(size)
            order.append(chrom)
    return sizes, order

def parse_attr(attr):
    d = {}
    for item in attr.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        m = re.match(r'(\S+)\s+"?([^"]+)"?', item)
        if m:
            d[m.group(1)] = m.group(2)
    return d

def safe_name(gene_id, gene_name):
    if gene_name and gene_name != gene_id:
        return f"{gene_id}|{gene_name}"
    return gene_id

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--promoter-up", type=int, default=1500)
    ap.add_argument("--promoter-down", type=int, default=1500)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    chrom_sizes, chrom_order = load_chrom_sizes(args.chrom_sizes)
    keep_chroms = set(chrom_sizes)

    genes = []

    with open_text(args.gtf) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attr = parts

            if feature != "gene":
                continue
            if chrom not in keep_chroms:
                continue
            if strand not in ["+", "-"]:
                continue

            start = int(start)
            end = int(end)

            # GTF: 1-based closed
            # BED: 0-based half-open
            bed_start = start - 1
            bed_end = end

            attrs = parse_attr(attr)
            gene_id = attrs.get("gene_id", f"{chrom}:{start}-{end}")
            gene_name = attrs.get("gene_name", gene_id)
            gene_biotype = (
                attrs.get("gene_biotype")
                or attrs.get("gene_type")
                or attrs.get("biotype")
                or "NA"
            )

            if bed_end <= bed_start:
                continue

            genes.append({
                "chrom": chrom,
                "start": bed_start,
                "end": bed_end,
                "strand": strand,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "gene_biotype": gene_biotype,
                "name": safe_name(gene_id, gene_name)
            })

    def chrom_key(g):
        try:
            cidx = chrom_order.index(g["chrom"])
        except ValueError:
            cidx = 999999
        return (cidx, g["start"], g["end"], g["gene_id"])

    genes.sort(key=chrom_key)

    gene_body = outdir / "pig_gene_body.bed"
    promoter = outdir / "pig_promoter_TSS_1.5kb.bed"
    gene_body_no_promoter = outdir / "pig_gene_body_no_promoter.bed"
    summary = outdir / "pig_gene_regions_summary.txt"

    n_gene = 0
    n_promoter = 0
    n_no_promoter = 0

    with open(gene_body, "w") as gb, open(promoter, "w") as pr, open(gene_body_no_promoter, "w") as gbnp:
        for g in genes:
            chrom = g["chrom"]
            start = g["start"]
            end = g["end"]
            strand = g["strand"]
            chrom_size = chrom_sizes[chrom]
            name = g["name"]

            length = end - start

            gb.write(
                f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\t"
                f"{g['gene_id']}\t{g['gene_name']}\t{g['gene_biotype']}\t{length}\n"
            )
            n_gene += 1

            # TSS coordinate in BED system
            # + strand: GTF start -> BED start
            # - strand: GTF end -> BED end coordinate
            if strand == "+":
                tss = start
            else:
                tss = end

            p_start = max(0, tss - args.promoter_up)
            p_end = min(chrom_size, tss + args.promoter_down)

            if p_end > p_start:
                pr.write(
                    f"{chrom}\t{p_start}\t{p_end}\t{name}\t0\t{strand}\t"
                    f"{g['gene_id']}\t{g['gene_name']}\t{g['gene_biotype']}\t{tss}\n"
                )
                n_promoter += 1

            # gene body without promoter-overlapping part
            if strand == "+":
                np_start = min(end, start + args.promoter_down)
                np_end = end
            else:
                np_start = start
                np_end = max(start, end - args.promoter_up)

            if np_end > np_start:
                np_len = np_end - np_start
                gbnp.write(
                    f"{chrom}\t{np_start}\t{np_end}\t{name}\t0\t{strand}\t"
                    f"{g['gene_id']}\t{g['gene_name']}\t{g['gene_biotype']}\t{np_len}\n"
                )
                n_no_promoter += 1

    with open(summary, "w") as fo:
        fo.write(f"gtf\t{args.gtf}\n")
        fo.write(f"chrom_sizes\t{args.chrom_sizes}\n")
        fo.write(f"promoter_up\t{args.promoter_up}\n")
        fo.write(f"promoter_down\t{args.promoter_down}\n")
        fo.write(f"genes_total\t{n_gene}\n")
        fo.write(f"promoters_total\t{n_promoter}\n")
        fo.write(f"gene_body_no_promoter_total\t{n_no_promoter}\n")

    print(f"[INFO] wrote {gene_body}")
    print(f"[INFO] wrote {promoter}")
    print(f"[INFO] wrote {gene_body_no_promoter}")
    print(f"[INFO] wrote {summary}")

if __name__ == "__main__":
    main()
