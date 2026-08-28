#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import plotly.express as px

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umap", required=True, help="UMAP TSV with UMAP1, UMAP2, UMAP3")
    ap.add_argument("--out", required=True, help="Output interactive html")
    ap.add_argument("--color", default="cluster", help="Column used for color")
    ap.add_argument("--title", default="3D promoter residual UMAP")
    ap.add_argument("--max-cells", type=int, default=80000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--marker-size", type=float, default=2.0)
    args = ap.parse_args()

    umap_file = Path(args.umap)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(umap_file, sep="\t")

    required = ["UMAP1", "UMAP2", "UMAP3"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"{c} not found in {umap_file}")

    if args.color not in df.columns:
        raise ValueError(f"{args.color} not found in {umap_file}. Available columns: {list(df.columns)}")

    # 交互网页细胞太多会卡，超过上限则只为网页展示随机抽样；原始tsv不变
    if len(df) > args.max_cells:
        df = df.sample(args.max_cells, random_state=args.seed).copy()

    # cluster作为离散分类显示；QC指标作为连续变量显示
    if args.color == "cluster":
        df[args.color] = df[args.color].astype(str)

    hover_cols = []
    for c in ["cluster", "reads_used", "covered_calls", "covered_bins", "modified_rate"]:
        if c in df.columns and c != args.color:
            hover_cols.append(c)

    fig = px.scatter_3d(
        df,
        x="UMAP1",
        y="UMAP2",
        z="UMAP3",
        color=args.color,
        hover_data=hover_cols,
        title=args.title,
        opacity=0.75,
    )

    fig.update_traces(marker=dict(size=args.marker_size))
    fig.update_layout(
        template="simple_white",
        scene=dict(
            xaxis_title="UMAP1",
            yaxis_title="UMAP2",
            zaxis_title="UMAP3",
        ),
        legend_title_text=args.color,
        margin=dict(l=0, r=0, b=0, t=50),
    )

    # inline模式生成自包含HTML，拷到本地也能打开旋转
    fig.write_html(out, include_plotlyjs=True, full_html=True)
    print("saved:", out)

if __name__ == "__main__":
    main()
