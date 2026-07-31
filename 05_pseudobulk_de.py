#!/usr/bin/env python
"""
05_pseudobulk_de.py

This script:
1. For each cell type (Memory B or Plasma) finds n = the minimum number of cells
   any single donor has for that cell type.
2. Randomly samples exactly n cells from each donor for that cell type
3. Sums raw counts per donor x cell type into one pseudobulk profile per
   donor/cell type (16 profiles total: 8 donors x 2 cell types).
4. Run pseudobulk differential expression (PyDESeq2) Memory_B vs Plasma,
   paired by donor.

Usage:
    python 05_pseudobulk_de.py --combined_h5ad ./celltype_out/combined_annotated.h5ad \
                                 --outdir ./pseudobulk_out --seed 0
"""
import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def equal_n_sample(adata, seed=0):
    rng = np.random.default_rng(seed)
    sampled_obs = []
    for cell_type in ["Memory_B", "Plasma"]:
        sub = adata.obs[adata.obs.cell_type == cell_type]
        n_per_donor = sub.groupby("donor").size()
        n = int(n_per_donor.min())   # the binding constraint
        print(f"{cell_type}: min donor count = {n} (from donor {n_per_donor.idxmin()}); "
              f"sampling {n} cells from every donor")
        for donor, grp in sub.groupby("donor"):
            chosen = rng.choice(grp.index.values, size=n, replace=False)
            sampled_obs.extend(chosen)
    return adata[adata.obs_names.isin(sampled_obs)].copy()


def make_pseudobulk(adata):
    rows = []
    meta = []
    for (donor, cell_type), grp in adata.obs.groupby(["donor", "cell_type"]):
        if len(grp) == 0:
            continue
        cell_idx = grp.index
        counts = adata[cell_idx].layers["counts"]
        summed = np.asarray(counts.sum(axis=0)).ravel()
        rows.append(summed)
        meta.append({"donor": donor, "cell_type": cell_type,
                      "sex": grp["sex"].iloc[0], "age": grp["age"].iloc[0],
                      "n_cells": len(grp)})
    counts_df = pd.DataFrame(rows, columns=adata.var_names,
                              index=[f"{m['donor']}_{m['cell_type']}" for m in meta])
    meta_df = pd.DataFrame(meta, index=counts_df.index)
    return counts_df, meta_df


def main(combined_h5ad, outdir, seed):
    os.makedirs(outdir, exist_ok=True)
    adata = sc.read_h5ad(combined_h5ad)
    adata = adata[adata.obs.cell_type.isin(["Memory_B", "Plasma"])].copy()

    balanced = equal_n_sample(adata, seed=seed)
    balanced.write(f"{outdir}/balanced_cells.h5ad")

    counts_df, meta_df = make_pseudobulk(balanced)
    counts_df = counts_df.loc[:, counts_df.sum(axis=0) > 0]  # remove all-zero genes
    counts_df.to_csv(f"{outdir}/pseudobulk_counts.csv")
    meta_df.to_csv(f"{outdir}/pseudobulk_meta.csv")

    #DE: Memory_B vs Plasma, paired by donor
    meta_df["cell_type"] = pd.Categorical(meta_df["cell_type"],
                                           categories=["Plasma", "Memory_B"])  # Plasma as reference
    dds = DeseqDataSet(counts=counts_df.astype(int), metadata=meta_df,
                        design_factors=["donor", "cell_type"])
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["cell_type", "Memory_B", "Plasma"])
    stats.summary()
    res = stats.results_df.sort_values("padj")
    res.to_csv(f"{outdir}/de_memoryB_vs_plasma.csv")
    print(res.head(20))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--combined_h5ad", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    main(args.combined_h5ad, args.outdir, args.seed)
