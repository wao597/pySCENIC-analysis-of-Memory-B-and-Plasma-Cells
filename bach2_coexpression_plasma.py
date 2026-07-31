#!/usr/bin/env python
"""
bach2_coexpression_plasma.py

This script produces a data-derived figure showing BACH2 co-expression with TP63, RORA, FGF13
and FUT8 within plasma cells from the balanced_cells.h5ad. This directly visualises the Lasso inference; cells with higher BACH2 co-express TP63, RORA, FGF13 (positive predictors) and have lower FUT8
(negative predictor / glycosylation).

Usage:
    python bach2_coexpression_plasma.py \
        --h5ad /path/to/pseudobulk_out/balanced_cells.h5ad \
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
import matplotlib.gridspec as gridspec
from scipy.stats import pearsonr, spearmanr

BLUE   = "#1565C0"
RED    = "#C62828"
ORANGE = "#E65100"
GREY   = "#607D8B"

# genes to plot against BACH2
# positive Lasso predictors in plasma cells
POSITIVE_GENES = ["TP63", "RORA", "FGF13"]
# negative Lasso predictor 
NEGATIVE_GENE  = "FUT8"


def get_expr(adata, gene):
    """Extract log-normalised expression for a gene as a numpy array."""
    if gene not in adata.var_names:
        return None
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x).ravel()


def scatter_with_trend(ax, x, y, gene_x, gene_y, colour, title):
    """Scatter plot with regression line and Spearman r.
    Only computes correlation on cells expressing BACH2 (x > 0)
    to avoid zero-inflation bias from double-zero cells."""

    # all cells
    ax.scatter(x, y, c="#ECEFF1", alpha=0.3, s=4,
               linewidths=0, rasterized=True, zorder=1)

    # cells with BACH2 expression > 0 — coloured
    mask = x > 0
    x_m, y_m = x[mask], y[mask]

    # colour by BACH2 level
    if len(x_m) > 0:
        ax.scatter(x_m, y_m, c=x_m, cmap="YlOrRd",
                   alpha=0.6, s=10, linewidths=0,
                   rasterized=True, zorder=2,
                   vmin=x_m.min(), vmax=x_m.max())

    # regression line on BACH2-expressing cells only
    if len(x_m) > 10:
        m, b = np.polyfit(x_m, y_m, 1)
        xline = np.linspace(x_m.min(), x_m.max(), 100)
        ax.plot(xline, m * xline + b, color="black", lw=1.8, zorder=3)

    # Spearman correlation on BACH2-expressing cells
    if len(x_m) > 10:
        rho, pval = spearmanr(x_m, y_m)
        pstr = f"p={pval:.2e}" if pval >= 1e-4 else "p<1e-4"
        col = "#C62828" if rho < 0 else "#1B5E20"
        ax.text(0.05, 0.93, f"ρ={rho:.2f}, {pstr}\n(BACH2+ cells, n={len(x_m)})",
                transform=ax.transAxes, fontsize=8,
                color=col, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white",
                          ec="lightgrey", alpha=0.9))

    ax.set_xlabel("BACH2 expression (log-normalised)", fontsize=9)
    ax.set_ylabel(f"{gene_y} expression (log-normalised)", fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold", pad=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main(h5ad, outdir):
    os.makedirs(outdir, exist_ok=True)

    print("Loading balanced cells...", flush=True)
    adata = sc.read_h5ad(h5ad)

    # subset to plasma cells
    plasma = adata[adata.obs["cell_type"] == "Plasma"].copy()
    print(f"Plasma cells: {plasma.n_obs}", flush=True)

    # get BACH2 expression
    bach2 = get_expr(plasma, "BACH2")
    if bach2 is None:
        raise ValueError("BACH2 not found in expression matrix")

    # check which genes are available
    genes_pos = [g for g in POSITIVE_GENES if g in plasma.var_names]
    gene_neg  = NEGATIVE_GENE if NEGATIVE_GENE in plasma.var_names else None

    missing = [g for g in POSITIVE_GENES if g not in plasma.var_names]
    if missing:
        print(f"WARNING: genes not found, will skip: {missing}", flush=True)

    n_panels = len(genes_pos) + (1 if gene_neg else 0)
    if n_panels == 0:
        raise ValueError("No target genes found in expression matrix")

    
    n_cols = 2
    n_rows = int(np.ceil(n_panels / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(11, 5 * n_rows),
                              facecolor="white")
    axes = np.atleast_1d(axes).ravel()

    # positive predictor panels
    for i, gene in enumerate(genes_pos):
        y = get_expr(plasma, gene)
        scatter_with_trend(axes[i], bach2, y, "BACH2", gene,
                           colour=ORANGE,
                           title=f"BACH2 vs {gene}\n(positive predictor in plasma)")

    # negative predictor panel
    if gene_neg:
        y = get_expr(plasma, gene_neg)
        scatter_with_trend(axes[len(genes_pos)], bach2, y, "BACH2", gene_neg,
                           colour=RED,
                           title=f"BACH2 vs {gene_neg}\n(negative predictor — glycosylation)")

    #contigency for unused axes
    for j in range(n_panels, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "BACH2 co-expression within Plasma cells\n"
        "Cells with higher BACH2 co-express TP63, RORA & FGF13 (positive)\n"
        "and show lower FUT8/glycosylation activity (negative)",
        fontsize=11, fontweight="bold", y=1.02
    )
    fig.tight_layout()

    out = os.path.join(outdir, "bach2_coexpression_plasma.png")
    fig.savefig(out, dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out}", flush=True)

    # Expression statistics
    print("\nExpression statistics in plasma cells:")
    for gene in genes_pos + ([gene_neg] if gene_neg else []):
        y = get_expr(plasma, gene)
        n_expr = (y > 0).sum()
        print(f"  {gene}: {n_expr}/{plasma.n_obs} cells expressing "
              f"({100*n_expr/plasma.n_obs:.1f}%)", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad",   required=True,
                   help="Path to balanced_cells.h5ad from pseudobulk step")
    p.add_argument("--outdir", required=True)
    args = p.parse_args()
    main(args.h5ad, args.outdir)
