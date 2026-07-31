#!/usr/bin/env python
"""
top_de_barplot.py

This scripts produces a horizontal bar chart of top 20 DE genes per cell type.


Usage:
    python top_de_barplot.py \
        --de_csv /path/to/de_memoryB_vs_plasma.csv \
        --outdir /path/to/output \
        --top_n 10
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def make_barplot(de_csv, outdir, top_n=10):
    os.makedirs(outdir, exist_ok=True)

    de = pd.read_csv(de_csv, index_col=0)
    de = de.dropna(subset=["padj", "log2FoldChange"])
    sig = de[(de["padj"] < 0.05) & (de["log2FoldChange"].abs() > 1)]

    # top N per direction by absolute log2FC
    top_mb = sig[sig["log2FoldChange"] > 0].nlargest(top_n, "log2FoldChange")
    top_pl = sig[sig["log2FoldChange"] < 0].nsmallest(top_n, "log2FoldChange")

    # combine plasma on top (negative), memory B below (positive) sort so most extreme are nearest the centre
    top_pl_sorted = top_pl.sort_values("log2FoldChange", ascending=False)
    top_mb_sorted = top_mb.sort_values("log2FoldChange", ascending=True)
    combined = pd.concat([top_pl_sorted, top_mb_sorted])

    colours = ["#C62828" if v < 0 else "#1565C0"
               for v in combined["log2FoldChange"]]

    fig, ax = plt.subplots(figsize=(9, 8), facecolor="white")

    bars = ax.barh(combined.index, combined["log2FoldChange"],
                   color=colours, edgecolor="white", height=0.7, alpha=0.88)

    # add padj stars
    for i, (gene, row) in enumerate(combined.iterrows()):
        x = row["log2FoldChange"]
        padj = row["padj"]
        stars = "***" if padj < 1e-10 else "**" if padj < 1e-5 else "*"
        x_pos = x + (0.15 if x > 0 else -0.15)
        ha = "left" if x > 0 else "right"
        ax.text(x_pos, i, stars, va="center", ha=ha,
                fontsize=8, color="#37474F", fontweight="bold")

    ax.axvline(0, color="black", linewidth=1.0, zorder=3)

    # shade regions
    xlim = ax.get_xlim()
    ax.axvspan(xlim[0], 0, alpha=0.04, color="#C62828", zorder=0)
    ax.axvspan(0, xlim[1], alpha=0.04, color="#1565C0", zorder=0)

    # cell type labels
    ax.text(xlim[0] + 0.1, top_n + 0.2, "← Higher in Plasma",
            fontsize=10, color="#C62828", fontweight="bold", fontstyle="italic")
    ax.text(xlim[1] - 0.1, top_n + 0.2, "Higher in Memory B →",
            fontsize=10, color="#1565C0", fontweight="bold",
            fontstyle="italic", ha="right")

    ax.set_xlabel("log₂ Fold Change (Memory B vs Plasma)",
                  fontsize=12, fontweight="bold")
    ax.set_title(
        f"Top {top_n} Differentially Expressed Genes per Cell Type\n"
        "Memory B cells vs Plasma cells · PyDESeq2 · padj < 0.05 · |log₂FC| > 1",
        fontsize=11, fontweight="bold", pad=12)

    ax.tick_params(axis="y", labelsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.yaxis.set_tick_params(length=0)

    # colour y-tick labels to match bar colour
    for tick, gene in zip(ax.get_yticklabels(), combined.index):
        fc = combined.loc[gene, "log2FoldChange"]
        tick.set_color("#C62828" if fc < 0 else "#1565C0")
        tick.set_fontweight("bold")

    fig.tight_layout()
    out = f"{outdir}/top_de_barplot.png"
    fig.savefig(out, dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out}")
    print("\nGenes plotted:")
    print(combined[["log2FoldChange", "padj"]].to_string())

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--de_csv",  required=True)
    p.add_argument("--outdir",  required=True)
    p.add_argument("--top_n",   type=int, default=10)
    args = p.parse_args()
    make_barplot(args.de_csv, args.outdir, args.top_n)
