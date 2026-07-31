#!/usr/bin/env python
"""
bach2_high_low_violin.py

This script divides plasma cells into BACH2-high and BACH2-low groups (top/bottom tertile of BACH2-expressing cells) and compares expression of TP63, RORA, FGF13 (positive Lasso predictors) and FUT8
(negative predictor) between groups using violin plots with
Mann-Whitney U test. This directly visualises whether BACH2-high plasma cells have
higher TP63/RORA/FGF13 and lower FUT8 than BACH2-low plasma cells.

Usage:
    python bach2_high_low_violin.py \
        --h5ad /path/to/balanced_cells.h5ad \
        --outdir /path/to/output
"""
import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ORANGE  = "#E65100"
BLUE    = "#1565C0"
RED     = "#C62828"
LGREY   = "#ECEFF1"

GENES = ["TP63", "RORA", "FGF13", "FUT8"]
COLOURS = {
    "BACH2-high": ORANGE,
    "BACH2-low":  BLUE,
}


def get_expr(adata, gene):
    if gene not in adata.var_names:
        return None
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x).ravel()


def stars(pval):
    if pval < 0.001: return "***"
    elif pval < 0.01: return "**"
    elif pval < 0.05: return "*"
    else: return "ns"


def main(h5ad, outdir):
    os.makedirs(outdir, exist_ok=True)

    print("Loading balanced cells...", flush=True)
    adata = sc.read_h5ad(h5ad)
    plasma = adata[adata.obs["cell_type"] == "Plasma"].copy()
    print(f"Total plasma cells: {plasma.n_obs}", flush=True)

    # get BACH2 expression
    bach2 = get_expr(plasma, "BACH2")

    # keep only BACH2-expressing cells for the grouping
    bach2_pos_mask = bach2 > 0
    plasma_pos = plasma[bach2_pos_mask].copy()
    bach2_pos = bach2[bach2_pos_mask]
    print(f"BACH2+ plasma cells: {plasma_pos.n_obs}", flush=True)

    # split into top and bottom tertile by BACH2 expression
    t33 = np.percentile(bach2_pos, 33)
    t67 = np.percentile(bach2_pos, 67)

    high_mask = bach2_pos >= t67
    low_mask  = bach2_pos <= t33

    plasma_high = plasma_pos[high_mask]
    plasma_low  = plasma_pos[low_mask]

    print(f"BACH2-high (top tertile, BACH2 >= {t67:.2f}): {plasma_high.n_obs} cells",
          flush=True)
    print(f"BACH2-low (bottom tertile, BACH2 <= {t33:.2f}): {plasma_low.n_obs} cells",
          flush=True)

    # check target genes available
    available = [g for g in GENES if g in plasma.var_names]
    missing   = [g for g in GENES if g not in plasma.var_names]
    if missing:
        print(f"WARNING: not found, skipping: {missing}", flush=True)

    n = len(available)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 5.5), facecolor="white")
    if n == 1:
        axes = [axes]

    for ax, gene in zip(axes, available):
        y_high = get_expr(plasma_high, gene)
        y_low  = get_expr(plasma_low,  gene)

        # Mann-Whitney U test
        stat, pval = mannwhitneyu(y_high, y_low, alternative="two-sided")
        sig = stars(pval)

        # violin plot
        parts = ax.violinplot([y_low, y_high],
                              positions=[1, 2],
                              showmedians=True,
                              showextrema=False)

        for pc, col in zip(parts["bodies"],
                           [COLOURS["BACH2-low"], COLOURS["BACH2-high"]]):
            pc.set_facecolor(col)
            pc.set_alpha(0.7)
            pc.set_edgecolor("white")
        parts["cmedians"].set_colors(["white", "white"])
        parts["cmedians"].set_linewidth(2)

        # individual points
        for pos, y, col in [(1, y_low, COLOURS["BACH2-low"]),
                            (2, y_high, COLOURS["BACH2-high"])]:
            jitter = np.random.default_rng(42).uniform(-0.06, 0.06, len(y))
            ax.scatter(pos + jitter, y, c=col, alpha=0.4, s=6,
                       linewidths=0, zorder=3)

        # significance brackets
        ymax = max(y_high.max(), y_low.max())
        bracket_y = ymax * 1.08
        ax.plot([1, 1, 2, 2],
                [ymax * 1.02, bracket_y, bracket_y, ymax * 1.02],
                c="black", lw=1.0)
        ax.text(1.5, bracket_y * 1.02, sig,
                ha="center", va="bottom", fontsize=13, fontweight="bold",
                color="black" if sig != "ns" else LGREY)

        # labels
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["BACH2-low\n(bottom tertile)",
                             "BACH2-high\n(top tertile)"],
                           fontsize=9)
        ax.set_ylabel(f"{gene} expression\n(log-normalised)", fontsize=9)
        ax.set_title(f"{gene}", fontsize=12, fontweight="bold", pad=6)

        # colour x tick labels
        for tick, col in zip(ax.get_xticklabels(),
                             [COLOURS["BACH2-low"], COLOURS["BACH2-high"]]):
            tick.set_color(col)
            tick.set_fontweight("bold")

        ax.text(0.97, 0.97,
                f"p={pval:.3f}" if pval >= 0.001 else f"p={pval:.2e}",
                transform=ax.transAxes, fontsize=8,
                color="black", ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.2", fc="white",
                          ec="lightgrey", alpha=0.9))

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(0.4, 2.6)

    fig.suptitle(
        "BACH2-high vs BACH2-low Plasma Cells\n"
        "Gene expression of Lasso-identified co-expression partners\n"
        f"(BACH2+ plasma cells only, n={plasma_pos.n_obs}; "
        f"top/bottom tertile split)",
        fontsize=11, fontweight="bold", y=1.02
    )
    fig.tight_layout()

    out = os.path.join(outdir, "bach2_high_low_violin.png")
    fig.savefig(out, dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\nSaved: {out}", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad",   required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()
    main(args.h5ad, args.outdir)
