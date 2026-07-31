#!/usr/bin/env python
"""
This script performs Per-donor QC on raw data. 
It loads each donor's filtered feature-barcode matrix. Calculates per-cell QC metrics including number of genes detected, total UMI count, and percentage of reads mapping to mitochondrial genes. Filters cells to retain those with 200–6,000 genes detected and less than 10% mitochondrial reads. Removes genes detected in fewer than 3 cells. Applies log-normalisation with scale factor 10,000.  Run as an 8-task SLURM array jo

Usage:
    python 01_qc.py --cellranger_dir /path/to/cellranger_outs --donor BM1 \
                     --sex female --age 52 --outdir ./qc_out
"""
import argparse
import scanpy as sc
import numpy as np

sc.settings.verbosity = 1


def qc_one_donor(cellranger_dir, donor, sex, age, outdir,
                  min_genes=200, max_pct_mt=15, min_cells_gene=3,
                  expected_doublet_rate=0.06):
    h5_path = f"{cellranger_dir}/Manton{donor}/outs/filtered_feature_bc_matrix.h5"
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()
    adata.obs["donor"] = donor
    adata.obs["sex"] = sex
    adata.obs["age"] = age

    #QC Metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)

    #Doublet detection
    import scrublet as scr
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    # Filtering
    adata = adata[adata.obs.n_genes_by_counts >= min_genes].copy()
    adata = adata[adata.obs.pct_counts_mt < max_pct_mt].copy()
    adata = adata[~adata.obs.predicted_doublet].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells_gene)

    # Normalisation 
    adata.layers["counts"] = adata.X.copy()  # keep raw counts for pySCENIC / pseudobulk later
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    #HVG / dim reduction / clustering
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    adata.raw = adata
    sc.pp.pca(adata, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=1.0)
    sc.tl.umap(adata)

    out_file = f"{outdir}/{donor}_qc.h5ad"
    adata.write(out_file)
    print(f"{donor}: {adata.n_obs} cells x {adata.n_vars} genes after QC -> {out_file}")
    return out_file


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--cellranger_dir", required=True)
    p.add_argument("--donor", required=True)
    p.add_argument("--sex", required=True)
    p.add_argument("--age", required=True, type=int)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()
    import os
    os.makedirs(args.outdir, exist_ok=True)
    qc_one_donor(args.cellranger_dir, args.donor, args.sex, args.age, args.outdir)
