#!/usr/bin/env python3
import argparse
import gzip
import csv
import subprocess
from pathlib import Path


def norm(x):
    return x.strip().lower().replace("-", "_").replace(".", "_")


def pick_col(header, candidates):
    hmap = {norm(h): h for h in header}
    for c in candidates:
        if c in hmap:
            return hmap[c]
    return None


def to_float(x):
    try:
        if x is None or x == "" or x.lower() == "nan":
            return None
        return float(x)
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--genome-size", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--min-cov", type=int, default=5)
    ap.add_argument("--bedGraphToBigWig", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw_bg = outdir / f"{args.prefix}.raw.bedGraph"
    sorted_bg = outdir / f"{args.prefix}.sorted.bedGraph"
    bw = outdir / f"{args.prefix}.bw"

    with gzip.open(args.input, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = reader.fieldnames
        if header is None:
            raise RuntimeError("Empty input file or no header found")

        print("[INFO] input columns:", header, flush=True)

        chrom_col = pick_col(header, ["chr", "chrom", "chromosome", "seqname", "seqnames"])
        pos_col = pick_col(header, ["pos", "position", "site", "cpg_pos", "start"])
        end_col = pick_col(header, ["end"])

        cov_col = pick_col(header, [
            "covered_calls", "total_covered_calls", "coverage", "cov",
            "n_cov", "total", "depth", "n_total"
        ])

        rate_col = pick_col(header, [
            "modified_rate", "methylation_rate", "meth_rate",
            "conversion_rate", "cpg_conversion_rate", "rate", "ratio"
        ])

        mod_col = pick_col(header, [
            "modified_calls", "total_modified_calls", "methylated_calls",
            "n_mod", "n_methylated", "converted_calls"
        ])

        if pos_col is None:
            raise RuntimeError(f"Cannot find position column in {header}")

        if rate_col is None and not (mod_col is not None and cov_col is not None):
            raise RuntimeError(
                "Cannot find methylation/conversion value column. "
                "Need rate column, or modified_calls + covered_calls."
            )

        print("[INFO] chrom_col =", chrom_col, flush=True)
        print("[INFO] pos_col   =", pos_col, flush=True)
        print("[INFO] cov_col   =", cov_col, flush=True)
        print("[INFO] rate_col  =", rate_col, flush=True)
        print("[INFO] mod_col   =", mod_col, flush=True)

        n_in = 0
        n_out = 0

        with open(raw_bg, "w") as out:
            for row in reader:
                n_in += 1

                cov = None
                if cov_col is not None:
                    cov = to_float(row.get(cov_col))
                    if cov is None or cov < args.min_cov:
                        continue

                pos = to_float(row.get(pos_col))
                if pos is None:
                    continue

                if rate_col is not None:
                    val = to_float(row.get(rate_col))
                else:
                    mod = to_float(row.get(mod_col))
                    if mod is None or cov is None or cov == 0:
                        continue
                    val = mod / cov

                if val is None:
                    continue

                # 如果是0-1比例，转成0-100；如果本来就是百分比，保持原值
                if val <= 1.5:
                    val = val * 100.0

                chrom = row.get(chrom_col) if chrom_col is not None else args.chrom
                if chrom is None or chrom == "":
                    chrom = args.chrom

                # CpG位点通常是1-based单点；bedGraph需要0-based half-open
                start = int(pos) - 1
                end = int(pos)

                if start < 0 or end <= start:
                    continue

                out.write(f"{chrom}\t{start}\t{end}\t{val:.6f}\n")
                n_out += 1

    print("[INFO] total input rows:", n_in, flush=True)
    print("[INFO] output bedGraph rows:", n_out, flush=True)

    if n_out == 0:
        raise RuntimeError("No bedGraph rows generated. Check min_cov, columns, or input file.")

    print("[INFO] sorting bedGraph", flush=True)
    with open(sorted_bg, "w") as fout:
        subprocess.check_call(["sort", "-k1,1", "-k2,2n", str(raw_bg)], stdout=fout)

    print("[INFO] converting to bigWig", flush=True)
    subprocess.check_call([args.bedGraphToBigWig, str(sorted_bg), args.genome_size, str(bw)])

    print("[DONE]", bw, flush=True)


if __name__ == "__main__":
    main()
