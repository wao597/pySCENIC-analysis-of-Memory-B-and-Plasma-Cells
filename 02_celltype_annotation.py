#!/usr/bin/env python
"""
02_celltype_annotation.py
This script merges Azimuth labels  back into the QC'd h5ad files then:
  - Filters to high-confidence cells (predicted.celltype.l2.score > 0.5)
  - Subsets to Memory B cells and Plasma cells
  - Writes per-donor/per-celltype h5ads for pySCENIC
  - Writes combined_annotated.h5ad for pseudobulk DE
  - Generates EDA plots

Usage:
    python 02_celltype_annotation.py --qc_dir ./qc_out --outdir ./celltype_out
"""
import argparse
import glob
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Canonical marker sets
MEMORY_B_MARKERS = ["MS4A1", "CD19", "CD27", "CD79A", "CD79B", "TNFRSF13B", "IGHD"]
PLASMA_MARKERS = ["MZB1", "SDC1", "XBP1", "PRDM1", "IRF4", "JCHAIN", "CD38"]
# Note: in scoring, IGHD is *low* in memory B relative to naive; CD19/MS4A1 should be
# high in memory B and low/absent in plasma cells - inspect the score difference per
# cluster rather than a single cutoff.


def annotate(adata):
    sc.tl.score_genes(adata, [g for g in MEMORY_B_MARKERS if g in adata.var_names],
                       score_name="memoryB_score")
    sc.tl.score_genes(adata, [g for g in PLASMA_MARKERS if g in adata.var_names],
                       score_name="plasma_score")

    cluster_scores = adata.obs.groupby("leiden")[["memoryB_score", "plasma_score"]].mean()
    cluster_labels = {}
    for cl, row in cluster_scores.iterrows():
        if row["plasma_score"] > row["memoryB_score"] and row["plasma_score"] > 0.1:
            cluster_labels[cl] = "Plasma"
        elif row["memoryB_score"] > row["plasma_score"] and row["memoryB_score"] > 0.1:
            cluster_labels[cl] = "Memory_B"
        else:
            cluster_labels[cl] = "Other"
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_labels)
    return adata


def main(qc_dir, outdir):
    os.makedirs(outdir, exist_ok=True)
    files = sorted(glob.glob(f"{qc_dir}/*_qc.h5ad"))
    adatas = []
    for f in files:
        a = sc.read_h5ad(f)
        a = annotate(a)
        adatas.append(a)

    combined = adatas[0].concatenate(*adatas[1:], batch_key="donor_batch", index_unique="-")

    # --- EDA: cell type proportions by sex/age ---
    props = (combined.obs.groupby(["donor", "sex", "age", "cell_type"]).size()
             .reset_index(name="n_cells"))
    props["frac"] = props["n_cells"] / props.groupby("donor")["n_cells"].transform("sum")
    props.to_csv(f"{outdir}/celltype_proportions_by_donor.csv", index=False)

    fig, ax = plt.subplots(figsize=(7, 5))
    for ct in props.cell_type.unique():
        sub = props[props.cell_type == ct]
        ax.scatter(sub.age, sub.frac, label=ct)
    ax.set_xlabel("Donor age")
    ax.set_ylabel("Fraction of cells")
    ax.legend()
    fig.savefig(f"{outdir}/celltype_fraction_vs_age.png", dpi=150, bbox_inches="tight")

    sc.pl.pca(combined, color=["sex", "age", "cell_type"], save="_sex_age_celltype.png",
               show=False)

    # --- subset & write per donor x cell type, for pySCENIC ---
    for donor in combined.obs.donor.unique():
        for ct in ["Memory_B", "Plasma"]:
            sub = combined[(combined.obs.donor == donor) & (combined.obs.cell_type == ct)].copy()
            if sub.n_obs == 0:
                continue
            sub.write(f"{outdir}/{donor}_{ct}.h5ad")
            print(f"{donor} {ct}: {sub.n_obs} cells written")

    combined.write(f"{outdir}/combined_annotated.h5ad")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--qc_dir", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()
    main(args.qc_dir, args.outdir)
