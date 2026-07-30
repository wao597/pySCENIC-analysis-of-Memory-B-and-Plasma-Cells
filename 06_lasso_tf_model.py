#!/usr/bin/env python
"""
06_lasso_tf_model.py

For a chosen TF of interest (pick one of the top recurrent regulators from
step 4), train a Lasso model PER CELL TYPE on the balanced/matched-size
cell data from step 5:

    target  y = AUCell regulon activity score for the TF (continuous), or
                log-normalised expression of the TF itself if you'd rather
                model expression directly
    features X = log-normalised expression of all other genes

This tells you which genes are most predictive of (i.e. co-vary with /
plausibly regulated alongside) that TF's activity, within each cell type
separately - so you can compare its regulatory footprint in Memory B cells
vs Plasma cells.

Usage:
    python 06_lasso_tf_model.py --balanced_h5ad ./pseudobulk_out/balanced_cells.h5ad \
                                 --auc_dir ./pyscenic_out \
                                 --tf PRDM1 \
                                 --outdir ./lasso_out
"""
import argparse
import glob
import os
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LassoCV
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error
import matplotlib.pyplot as plt


def load_auc_scores(auc_dir, tf):
    """Stitch together this TF's AUCell score for every donor/cell_type run."""
    frames = []
    for auc_path in glob.glob(f"{auc_dir}/*/auc_mtx.csv"):
        run_name = os.path.basename(os.path.dirname(auc_path))
        df = pd.read_csv(auc_path, index_col=0)
        regulon_col = [c for c in df.columns if c.startswith(f"{tf}(")]
        if not regulon_col:
            continue
        s = df[regulon_col[0]].rename("tf_activity")
        s.index = [f"{idx}-{run_name}" for idx in s.index]  # match h5ad's obs_names suffixing if needed
        frames.append(s)
    return pd.concat(frames)


def fit_lasso_for_celltype(adata_ct, tf, outdir, cell_type):
    # y: use TF's own expression as the modelling target (robust fallback if
    # AUCell index matching against obs_names is fiddly across loom exports)
    if tf not in adata_ct.var_names:
        raise ValueError(f"{tf} not found in expression matrix for {cell_type}")

    y = np.asarray(adata_ct[:, tf].X.todense()).ravel() if hasattr(adata_ct.X, "todense") \
        else np.asarray(adata_ct[:, tf].X).ravel()
    X = adata_ct[:, adata_ct.var_names != tf].X
    X = X.toarray() if hasattr(X, "toarray") else X
    gene_names = adata_ct.var_names[adata_ct.var_names != tf]

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    model = LassoCV(cv=5, n_jobs=-1, max_iter=10000, random_state=0)
    model.fit(X_train_s, y_train)

    y_pred = model.predict(X_test_s)
    r2 = r2_score(y_test, y_pred)
    mse = mean_squared_error(y_test, y_pred)

    coefs = pd.Series(model.coef_, index=gene_names).sort_values(key=abs, ascending=False)
    nonzero = coefs[coefs != 0]
    nonzero.to_csv(f"{outdir}/{cell_type}_{tf}_lasso_coefs.csv")

    with open(f"{outdir}/{cell_type}_{tf}_metrics.txt", "w") as fh:
        fh.write(f"alpha (lambda): {model.alpha_}\n")
        fh.write(f"test R2: {r2:.4f}\n")
        fh.write(f"test MSE: {mse:.4f}\n")
        fh.write(f"n_nonzero_features: {len(nonzero)} / {X.shape[1]}\n")

    # plot top features
    top = nonzero.head(20)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.barh(top.index[::-1], top.values[::-1])
    ax.set_xlabel("Lasso coefficient")
    ax.set_title(f"{cell_type}: top predictors of {tf} expression")
    fig.tight_layout()
    fig.savefig(f"{outdir}/{cell_type}_{tf}_top_features.png", dpi=200)
    plt.close(fig)

    print(f"[{cell_type}] {tf}: R2={r2:.3f}, MSE={mse:.3f}, {len(nonzero)} nonzero features")
    return model, nonzero, r2, mse


def main(balanced_h5ad, tf, outdir):
    os.makedirs(outdir, exist_ok=True)
    adata = sc.read_h5ad(balanced_h5ad)

    results = {}
    for cell_type in ["Memory_B", "Plasma"]:
        sub = adata[adata.obs.cell_type == cell_type].copy()
        results[cell_type] = fit_lasso_for_celltype(sub, tf, outdir, cell_type)

    # quick side-by-side comparison of which genes matter in each cell type
    mb_genes = set(results["Memory_B"][1].index)
    pc_genes = set(results["Plasma"][1].index)
    print(f"Shared predictive genes across both cell types: {mb_genes & pc_genes}")
    print(f"Memory_B-specific: {mb_genes - pc_genes}")
    print(f"Plasma-specific: {pc_genes - mb_genes}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--balanced_h5ad", required=True)
    p.add_argument("--auc_dir", required=False, help="not required for the expression-target fallback")
    p.add_argument("--tf", required=True, help="TF of interest, e.g. a top hit from step 4")
    p.add_argument("--outdir", required=True)
    args = p.parse_args()
    main(args.balanced_h5ad, args.tf, args.outdir)
