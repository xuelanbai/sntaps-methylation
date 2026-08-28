#!/usr/bin/env python3
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

all_qc = "08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.global_per_cb.tsv"

groups = {
    "all_CBs": all_qc,
    "expected_top20000": "08_qc_summary/final_single_sample/cell_filter_expected_real_cells/snTAPS-FL-1.expected_realcell_top20000.global_per_cb.tsv",
    "expected_top10000": "08_qc_summary/final_single_sample/cell_filter_expected_real_cells/snTAPS-FL-1.expected_realcell_top10000.global_per_cb.tsv",
    "expected_top5000": "08_qc_summary/final_single_sample/cell_filter_expected_real_cells/snTAPS-FL-1.expected_realcell_top5000.global_per_cb.tsv",
    "expected_top2000": "08_qc_summary/final_single_sample/cell_filter_expected_real_cells/snTAPS-FL-1.expected_realcell_top2000.global_per_cb.tsv",
    "balanced_core10000": "08_qc_summary/final_single_sample/cell_filter_balanced_core/snTAPS-FL-1.balanced_core10000.global_per_cb.tsv",
    "balanced_core5000": "08_qc_summary/final_single_sample/cell_filter_balanced_core/snTAPS-FL-1.balanced_core5000.global_per_cb.tsv",
    "balanced_core2000": "08_qc_summary/final_single_sample/cell_filter_balanced_core/snTAPS-FL-1.balanced_core2000.global_per_cb.tsv",
}

outdir = Path("08_qc_summary/final_single_sample/cell_filter_qc_report")
outdir.mkdir(parents=True, exist_ok=True)

summary_rows = []

for name, path in groups.items():
    p = Path(path)
    if not p.exists():
        print("[WARN] missing", path)
        continue

    df = pd.read_csv(p, sep="\t")
    for c in ["reads_used", "modified_calls", "covered_calls", "modified_rate", "covered_bins"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    row = {"group": name, "n_CBs": df.shape[0]}
    for c in ["reads_used", "covered_calls", "covered_bins", "modified_rate"]:
        if c in df.columns:
            x = df[c].dropna()
            row[f"{c}_median"] = x.median()
            row[f"{c}_q10"] = x.quantile(0.10)
            row[f"{c}_q90"] = x.quantile(0.90)
    summary_rows.append(row)

    # 分布图
    for c in ["reads_used", "covered_calls", "covered_bins", "modified_rate"]:
        if c not in df.columns:
            continue
        x = df[c].dropna()
        if len(x) == 0:
            continue

        plt.figure(figsize=(7, 5))
        if c == "modified_rate":
            plt.hist(x, bins=80, density=True, alpha=0.8)
            plt.xlabel(c)
        else:
            plt.hist(np.log10(x + 1), bins=80, density=True, alpha=0.8)
            plt.xlabel(f"log10({c}+1)")
        plt.ylabel("Density")
        plt.title(f"{name} {c} distribution")
        plt.axvline(np.median(np.log10(x + 1)) if c != "modified_rate" else np.median(x),
                    linestyle="--", label="median")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / f"{name}.{c}.density.png", dpi=200)
        plt.close()

# CB rank曲线只对all_CBs画
df = pd.read_csv(all_qc, sep="\t")
for c in ["reads_used", "covered_calls", "covered_bins"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")
    x = df[c].dropna().sort_values(ascending=False).to_numpy()
    rank = np.arange(1, len(x) + 1)

    plt.figure(figsize=(7, 5))
    plt.plot(np.log10(rank), np.log10(x + 1), linewidth=1)
    plt.xlabel("log10(CB rank)")
    plt.ylabel(f"log10({c}+1)")
    plt.title(f"CB rank curve by {c}")
    for n in [2000, 5000, 10000, 20000, 85014]:
        if n < len(x):
            plt.axvline(np.log10(n), linestyle="--")
            plt.text(np.log10(n), plt.ylim()[1]*0.95, str(n), rotation=90, va="top")
    plt.tight_layout()
    plt.savefig(outdir / f"all_CBs.rank_curve.{c}.png", dpi=200)
    plt.close()

summary = pd.DataFrame(summary_rows)
summary.to_csv(outdir / "cell_filter_qc_summary.tsv", sep="\t", index=False)
print("[INFO] wrote", outdir / "cell_filter_qc_summary.tsv")
