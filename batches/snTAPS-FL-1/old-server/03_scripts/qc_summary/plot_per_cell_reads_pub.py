import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

inp = Path("08_qc_summary/final_single_sample/snTAPS-FL-1.all_technical_merged.global_per_cb.tsv")
out = Path("08_qc_summary/final_single_sample/qc_read_stats")
out.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(inp, sep="\t")
reads = pd.to_numeric(df["reads_used"], errors="coerce").dropna()
reads = reads[reads > 0].astype(float)

# summary
q25 = reads.quantile(0.25)
q50 = reads.quantile(0.50)
q75 = reads.quantile(0.75)
q05 = reads.quantile(0.05)
q95 = reads.quantile(0.95)

summary = pd.DataFrame({
    "metric": ["n_cells", "mean", "median", "q05", "q25", "q75", "q95", "min", "max"],
    "value": [
        len(reads),
        reads.mean(),
        q50,
        q05,
        q25,
        q75,
        q95,
        reads.min(),
        reads.max()
    ]
})
summary.to_csv(out / "per_cell_reads_used_pub_summary.tsv", sep="\t", index=False)

# -----------------------------
# 1. publication-style density plot
# -----------------------------
x = np.log10(reads + 1)

bins = 120
hist, edges = np.histogram(x, bins=bins, density=True)
centers = (edges[:-1] + edges[1:]) / 2

# gaussian-like smoothing without scipy
def smooth_curve(y, window=11, sigma=2.5):
    if window % 2 == 0:
        window += 1
    k = np.arange(window) - window // 2
    kernel = np.exp(-(k**2) / (2 * sigma**2))
    kernel = kernel / kernel.sum()
    return np.convolve(y, kernel, mode="same")

smooth = smooth_curve(hist, window=17, sigma=3.0)

x25 = np.log10(q25 + 1)
x50 = np.log10(q50 + 1)
x75 = np.log10(q75 + 1)

plt.figure(figsize=(6.2, 4.8))

# light histogram as background
plt.hist(
    x,
    bins=bins,
    density=True,
    alpha=0.28,
    edgecolor="none",
    label="Observed cells"
)

# smoothed density curve
plt.plot(
    centers,
    smooth,
    linewidth=2.2,
    label="Smoothed density"
)

# median and IQR
plt.axvline(x50, linestyle="--", linewidth=1.6, label=f"Median={int(round(q50)):,}")
plt.axvline(x25, linestyle=":", linewidth=1.2, label=f"Q1={int(round(q25)):,}")
plt.axvline(x75, linestyle=":", linewidth=1.2, label=f"Q3={int(round(q75)):,}")

plt.xlabel("log10(reads per cell + 1)")
plt.ylabel("Density")
plt.title("Per-cell reads distribution")
plt.legend(frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig(out / "per_cell_reads_used_density_log10_pub.png", dpi=300)
plt.savefig(out / "per_cell_reads_used_density_log10_pub.pdf")
plt.close()

# -----------------------------
# 2. rank curve with smoothing
# -----------------------------
ranked = reads.sort_values(ascending=False).reset_index(drop=True)
rank = np.arange(1, len(ranked) + 1)

# smooth in log space, then back-transform
log_reads = np.log10(ranked.values + 1)
window = max(101, int(len(log_reads) * 0.002))
if window % 2 == 0:
    window += 1

kernel = np.ones(window) / window
smooth_log_reads = np.convolve(log_reads, kernel, mode="same")

# edge修正，避免开头结尾被卷积拉低
half = window // 2
smooth_log_reads[:half] = log_reads[:half]
smooth_log_reads[-half:] = log_reads[-half:]

plt.figure(figsize=(6.2, 4.8))
plt.plot(rank, ranked.values, linewidth=0.8, alpha=0.35, label="Observed")
plt.plot(rank, 10**smooth_log_reads - 1, linewidth=2.0, label="Smoothed")
plt.axhline(q50, linestyle="--", linewidth=1.4, label=f"Median={int(round(q50)):,}")

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Cell barcode rank")
plt.ylabel("Reads per cell")
plt.title("Per-cell reads rank curve")
plt.legend(frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig(out / "per_cell_reads_used_rank_curve_pub.png", dpi=300)
plt.savefig(out / "per_cell_reads_used_rank_curve_pub.pdf")
plt.close()

print("[DONE] outputs:")
print(out / "per_cell_reads_used_pub_summary.tsv")
print(out / "per_cell_reads_used_density_log10_pub.png")
print(out / "per_cell_reads_used_density_log10_pub.pdf")
print(out / "per_cell_reads_used_rank_curve_pub.png")
print(out / "per_cell_reads_used_rank_curve_pub.pdf")
