#!/usr/bin/env python
"""
06_lasso_tf_model.py

This script trains a Lasso model per cell type (Memory_B, Plasma) on balanced cell data,
predicting a TF's log-normalised expression from all other genes.
Pre-filters to top N highly variable genes before fitting to keep runtime
manageable (top 3000 HVGs, excluding the TF itself).

Usage:
    python 06_lasso_tf_model.py --balanced_h5ad /path/to/balanced_cells.h5ad \
                                 --tf BACH2 \
                                 --outdir /path/to/lasso_out/BACH2 \
                                 --n_hvg 3000
"""
import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LassoCV
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def fit_lasso_for_celltype(adata_ct, tf, outdir, cell_type, n_hvg):
    if tf not in adata_ct.var_names:
        raise ValueError(f"{tf} not found in expression matrix for {cell_type}")

    print(f"[{cell_type}] {adata_ct.n_obs} cells x {adata_ct.n_vars} genes before HVG filter", flush=True)

    # pre-filter to top N HVGs to keep Lasso runtime manageable 
    
    print(f"[{cell_type}] Selecting top {n_hvg} highly variable genes...", flush=True)
    sc.pp.highly_variable_genes(adata_ct, n_top_genes=n_hvg, flavor="seurat")
    hvg_names = adata_ct.var_names[adata_ct.var["highly_variable"]]

    # ensure TF is included
    all_keep = hvg_names.tolist()
    if tf not in all_keep:
        all_keep.append(tf)
    adata_sub = adata_ct[:, all_keep].copy()
    print(f"[{cell_type}] {adata_sub.n_vars} genes after HVG filter (including {tf})", flush=True)

    # y: TF log-normalised expression
    print(f"[{cell_type}] Extracting target variable ({tf})...", flush=True)
    y = np.asarray(adata_sub[:, tf].X.todense()).ravel() if hasattr(adata_sub.X, "todense") \
        else np.asarray(adata_sub[:, tf].X).ravel()

    # X: other genes
    mask = np.array(adata_sub.var_names) != tf
    X = adata_sub[:, mask].X
    X = X.toarray() if hasattr(X, "toarray") else np.array(X)
    gene_names = adata_sub.var_names[mask]
    print(f"[{cell_type}] Feature matrix: {X.shape[0]} cells x {X.shape[1]} features", flush=True)

    # train/test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)
    print(f"[{cell_type}] Train: {X_train.shape[0]} cells, Test: {X_test.shape[0]} cells", flush=True)

    # scale
    print(f"[{cell_type}] Scaling features...", flush=True)
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    # fit LassoCV
    print(f"[{cell_type}] Fitting LassoCV (5-fold CV, n_jobs=-1)...", flush=True)
    model = LassoCV(cv=5, n_jobs=-1, max_iter=10000, random_state=0)
    model.fit(X_train_s, y_train)
    print(f"[{cell_type}] LassoCV done. Alpha selected: {model.alpha_:.6f}", flush=True)

    # evaluate
    y_pred = model.predict(X_test_s)
    r2 = r2_score(y_test, y_pred)
    mse = mean_squared_error(y_test, y_pred)

    coefs = pd.Series(model.coef_, index=gene_names).sort_values(key=abs, ascending=False)
    nonzero = coefs[coefs != 0]

    # save coefficients
    nonzero.to_csv(f"{outdir}/{cell_type}_{tf}_lasso_coefs.csv")

    # save metrics
    with open(f"{outdir}/{cell_type}_{tf}_metrics.txt", "w") as fh:
        fh.write(f"TF: {tf}\n")
        fh.write(f"cell_type: {cell_type}\n")
        fh.write(f"n_hvg_used: {len(all_keep)}\n")
        fh.write(f"alpha (lambda): {model.alpha_}\n")
        fh.write(f"test R2: {r2:.4f}\n")
        fh.write(f"test MSE: {mse:.4f}\n")
        fh.write(f"n_nonzero_features: {len(nonzero)} / {X.shape[1]}\n")
        fh.write(f"train_cells: {X_train.shape[0]}\n")
        fh.write(f"test_cells: {X_test.shape[0]}\n")

    # plot top 20 features
    top = nonzero.head(20)
    colours = ["steelblue" if v > 0 else "salmon" for v in top.values[::-1]]
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.barh(top.index[::-1], top.values[::-1], color=colours)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Lasso coefficient")
    ax.set_title(f"{cell_type}: top 20 predictors of {tf} expression\nR²={r2:.3f}, MSE={mse:.3f}, alpha={model.alpha_:.4f}")
    fig.tight_layout()
    fig.savefig(f"{outdir}/{cell_type}_{tf}_top_features.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"[{cell_type}] {tf}: R2={r2:.3f}, MSE={mse:.3f}, {len(nonzero)} nonzero features", flush=True)
    return model, nonzero, r2, mse


def main(balanced_h5ad, tf, outdir, n_hvg):
    os.makedirs(outdir, exist_ok=True)
    print(f"Loading balanced h5ad: {balanced_h5ad}", flush=True)
    adata = sc.read_h5ad(balanced_h5ad)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes", flush=True)
    adata = adata[adata.obs.cell_type.isin(["Memory_B", "Plasma"])].copy()
    print(f"After cell type subset: {adata.n_obs} cells", flush=True)

    results = {}
    for cell_type in ["Memory_B", "Plasma"]:
        sub = adata[adata.obs.cell_type == cell_type].copy()
        print(f"\n--- {cell_type}: {sub.n_obs} cells ---", flush=True)
        results[cell_type] = fit_lasso_for_celltype(sub, tf, outdir, cell_type, n_hvg)

    # cross-cell-type comparison
    mb_genes = set(results["Memory_B"][1].index)
    pc_genes = set(results["Plasma"][1].index)
    shared = mb_genes & pc_genes

    comparison = pd.DataFrame({
        "Memory_B_coef": results["Memory_B"][1],
        "Plasma_coef": results["Plasma"][1]
    }).fillna(0)
    comparison.to_csv(f"{outdir}/{tf}_cross_celltype_coefs.csv")

    print(f"\n=== Cross cell type comparison for {tf} ===", flush=True)
    print(f"Shared predictive genes: {len(shared)}", flush=True)
    print(f"Memory_B-specific: {len(mb_genes - pc_genes)}", flush=True)
    print(f"Plasma-specific: {len(pc_genes - mb_genes)}", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--balanced_h5ad", required=True)
    p.add_argument("--tf", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--n_hvg", type=int, default=3000,
                   help="Number of highly variable genes to use as features (default 3000)")
    args = p.parse_args()
    main(args.balanced_h5ad, args.tf, args.outdir, args.n_hvg)
