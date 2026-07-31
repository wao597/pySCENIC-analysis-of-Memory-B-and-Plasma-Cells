#!/usr/bin/env python
"""
bach2_plasma_umap.py

This script produces a UMAP of plasma cells only, coloured by expression of BACH2/TP63/RORA/FGF13/FUT8. If BACH2-high cells spatially overlap with TP63/RORA/FGF13-high cells
on the UMAP, this provides visual evidence of the subpopulation.

Usage:
    python bach2_plasma_umap.py \
        --h5ad /path/to/balanced_cells.h5ad \
        --outdir /path/to/output
"""
import argparse
import os
import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENES = ["BACH2", "TP63", "RORA", "FGF13", "FUT8"]

# colour maps — BACH2 in orange, others in their own maps
CMAPS = {
    "BACH2": "YlOrRd",
    "TP63":  "Blues",
    "RORA":  "Purples",
    "FGF13": "Greens",
    "FUT8":  "Reds",
}


def main(h5ad, outdir):
    os.makedirs(outdir, exist_ok=True)

    print("Loading balanced cells...", flush=True)
    adata = sc.read_h5ad(h5ad)

    # subset to plasma cells only
    plasma = adata[adata.obs["cell_type"] == "Plasma"].copy()
    print(f"Plasma cells: {plasma.n_obs}", flush=True)

    # compute UMAP on plasma cells specifically
    print("Computing plasma-specific UMAP...", flush=True)
    sc.pp.highly_variable_genes(plasma, n_top_genes=2000, flavor="seurat")
    sc.pp.pca(plasma, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(plasma, n_neighbors=15)
    sc.tl.umap(plasma, random_state=0)

    # check which genes are available
    available = [g for g in GENES if g in plasma.var_names]
    missing   = [g for g in GENES if g not in plasma.var_names]
    if missing:
        print(f"WARNING: not found, will skip: {missing}", flush=True)

    n = len(available)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(5.5 * ncols, 5.0 * nrows),
                              facecolor="white")
    axes = axes.flatten() if n > 1 else [axes]

    umap1 = plasma.obsm["X_umap"][:, 0]
    umap2 = plasma.obsm["X_umap"][:, 1]

    for i, gene in enumerate(available):
        ax = axes[i]
        expr = plasma[:, gene].X
        if hasattr(expr, "toarray"):
            expr = expr.toarray()
        expr = np.asarray(expr).ravel()

        # plot non-expressing cells first (grey)
        zero_mask = expr == 0
        ax.scatter(umap1[zero_mask], umap2[zero_mask],
                   c="#ECEFF1", s=6, alpha=0.4,
                   linewidths=0, rasterized=True, zorder=1)

        # plot expressing cells on top, coloured by expression level
        pos_mask = expr > 0
        if pos_mask.sum() > 0:
            sc_plot = ax.scatter(
                umap1[pos_mask], umap2[pos_mask],
                c=expr[pos_mask],
                cmap=CMAPS.get(gene, "viridis"),
                s=14, alpha=0.85,
                linewidths=0, rasterized=True, zorder=2,
                vmin=expr[pos_mask].min(),
                vmax=np.percentile(expr[pos_mask], 95)
            )
            plt.colorbar(sc_plot, ax=ax, shrink=0.6, pad=0.02,
                         label="log-normalised expression")

        n_expr = pos_mask.sum()
        pct = 100 * n_expr / plasma.n_obs
        ax.set_title(f"{gene}\n({n_expr} cells expressing, {pct:.1f}%)",
                     fontsize=11, fontweight="bold", pad=6,
                     color="darkorange" if gene == "BACH2" else "black")
        ax.set_xlabel("UMAP 1", fontsize=8)
        ax.set_ylabel("UMAP 2", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # hide unused axes
    for j in range(len(available), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "UMAP of Plasma Cells — BACH2 and Co-expression Partners\n"
        "Grey = non-expressing cells  ·  Coloured = expressing cells\n"
        "Spatial overlap of BACH2 with TP63/RORA/FGF13 indicates "
        "co-expressing subpopulation",
        fontsize=11, fontweight="bold", y=1.01
    )
    fig.tight_layout()

    out = os.path.join(outdir, "bach2_plasma_umap.png")
    fig.savefig(out, dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out}", flush=True)

    # Save UMAP
    fig2, ax2 = plt.subplots(figsize=(7, 6), facecolor="white")
    ax2.scatter(umap1[zero_mask], umap2[zero_mask],
                c="#ECEFF1", s=8, alpha=0.4,
                linewidths=0, rasterized=True, zorder=1)
    bach2_expr = plasma[:, "BACH2"].X
    if hasattr(bach2_expr, "toarray"):
        bach2_expr = bach2_expr.toarray()
    bach2_expr = np.asarray(bach2_expr).ravel()
    pos = bach2_expr > 0
    sc2 = ax2.scatter(umap1[pos], umap2[pos],
                      c=bach2_expr[pos],
                      cmap="YlOrRd", s=20, alpha=0.9,
                      linewidths=0, rasterized=True, zorder=2)
    plt.colorbar(sc2, ax=ax2, label="BACH2 log-normalised expression")
    ax2.set_title(f"BACH2 expression in Plasma cells\n"
                  f"({pos.sum()} BACH2+ cells, {100*pos.sum()/plasma.n_obs:.1f}%)",
                  fontsize=12, fontweight="bold")
    ax2.set_xlabel("UMAP 1"); ax2.set_ylabel("UMAP 2")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    fig2.tight_layout()
    out2 = os.path.join(outdir, "bach2_plasma_umap_solo.png")
    fig2.savefig(out2, dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(fig2)
    print(f"Saved: {out2}", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad",   required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()
    main(args.h5ad, args.outdir)
