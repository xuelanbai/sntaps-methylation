#!/usr/bin/env python3
import argparse, gzip, re

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def parse_attrs(attr):
    d = {}
    for part in attr.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        m = re.match(r'(\S+)\s+"([^"]+)"', part)
        if m:
            d[m.group(1)] = m.group(2)
        else:
            xs = part.split("=", 1)
            if len(xs) == 2:
                d[xs[0]] = xs[1]
    return d

ap = argparse.ArgumentParser()
ap.add_argument("--gtf", required=True)
ap.add_argument("--markers", required=True)
ap.add_argument("--outprefix", required=True)
ap.add_argument("--promoter", type=int, default=1500)
ap.add_argument("--add-chr", action="store_true", help="Add chr prefix to chromosome names")
args = ap.parse_args()

markers = []
with open(args.markers) as f:
    for line in f:
        x = line.strip()
        if x:
            markers.append(x)

marker_set = {x.upper() for x in markers}

genes = {}
with open_maybe_gz(args.gtf) as f:
    for line in f:
        if line.startswith("#"):
            continue
        arr = line.rstrip("\n").split("\t")
        if len(arr) < 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attr = arr
        if feature != "gene":
            continue

        attrs = parse_attrs(attr)
        gene_name = attrs.get("gene_name") or attrs.get("gene") or attrs.get("Name") or attrs.get("gene_id")
        gene_id = attrs.get("gene_id", gene_name)
        if not gene_name:
            continue

        key = gene_name.upper()
        if key not in marker_set:
            continue

        if args.add_chr and not chrom.startswith("chr"):
            chrom = "chr" + chrom

        s = int(start) - 1   # BED 0-based
        e = int(end)         # BED end
        genes[key] = {
            "chrom": chrom,
            "start": s,
            "end": e,
            "name": gene_name,
            "gene_id": gene_id,
            "strand": strand
        }

with open(args.outprefix + ".gene_body.bed", "w") as body, \
     open(args.outprefix + ".promoter_1p5kb.bed", "w") as prom, \
     open(args.outprefix + ".tss_1bp.bed", "w") as tss, \
     open(args.outprefix + ".found.tsv", "w") as found, \
     open(args.outprefix + ".missing.txt", "w") as missing:

    found.write("marker\tchrom\tstart_1based\tend_1based\tstrand\tgene_id\n")

    for m in markers:
        key = m.upper()
        if key not in genes:
            missing.write(m + "\n")
            continue

        g = genes[key]
        chrom, s, e, name, strand, gene_id = g["chrom"], g["start"], g["end"], g["name"], g["strand"], g["gene_id"]

        body.write(f"{chrom}\t{s}\t{e}\t{name}|gene_body|{strand}\t0\t{strand}\n")

        if strand == "+":
            tss0 = s
        elif strand == "-":
            tss0 = e - 1
        else:
            tss0 = s

        ps = max(0, tss0 - args.promoter)
        pe = tss0 + args.promoter + 1

        prom.write(f"{chrom}\t{ps}\t{pe}\t{name}|TSS±{args.promoter}bp|{strand}\t0\t{strand}\n")
        tss.write(f"{chrom}\t{tss0}\t{tss0+1}\t{name}|TSS|{strand}\t0\t{strand}\n")
        found.write(f"{name}\t{chrom}\t{s+1}\t{e}\t{strand}\t{gene_id}\n")

print("[DONE]")
print("Found:", len(genes))
print("Output prefix:", args.outprefix)
