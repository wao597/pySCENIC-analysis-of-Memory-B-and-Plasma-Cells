#!/usr/bin/env python
"""
02b_merge_azimuth_labels.py
This script merges Azimuth labels  back into the QC'd h5ad files then:
  - Filters to high-confidence cells (predicted.celltype.l2.score > 0.5)
  - Subsets to Memory B cells and Plasma cells
  - Writes per-donor/per-celltype h5ads for pySCENIC
  - Writes combined_annotated.h5ad for pseudobulk DE
  - Generates EDA plots

Usage:
    python 02b_merge_azimuth_labels.py --qc_dir /path/to/qc_all_samples \
                                        --azimuth_dir /path/to/azimuth_out \
                                        --outdir /path/to/celltype_out \
                                        --min_score 0.5
"""
import argparse
import glob
import os
import scanpy as sc
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Azimuth l2 label 
LABEL_MAP = {
    "Memory B"  : "Memory_B",
    "Plasma"    : "Plasma",
    "Plasmablast": "Plasma",   # plasmablasts grouped with plasma cells
}

DONOR_META = {
    "BM1": {"sex": "female", "age": 52},
    "BM2": {"sex": "male",   "age": 50},
    "BM3": {"sex": "male",   "age": 39},
    "BM4": {"sex": "male",   "age": 29},
    "BM5": {"sex": "male",   "age": 29},
    "BM6": {"sex": "female", "age": 26},
    "BM7": {"sex": "female", "age": 36},
    "BM8": {"sex": "female", "age": 32},
}


def main(qc_dir, azimuth_dir, outdir, min_score):
    os.makedirs(outdir, exist_ok=True)
    sc.settings.figdir = outdir

    adatas = []

    for donor, meta in DONOR_META.items():
        qc_path  = f"{qc_dir}/{donor}_qc.h5ad"
        csv_path = f"{azimuth_dir}/{donor}_azimuth_labels.csv"

        if not os.path.exists(qc_path):
            print(f"WARNING: QC h5ad not found for {donor}, skipping")
            continue
        if not os.path.exists(csv_path):
            print(f"WARNING: Azimuth labels not found for {donor}, skipping")
            continue

        print(f"\n=== {donor} ===")
        adata = sc.read_h5ad(qc_path)
        labels = pd.read_csv(csv_path, index_col="cell_id")

        # align barcodes - QC h5ad barcodes may have donor suffix
        
        common = adata.obs_names.intersection(labels.index)
        if len(common) == 0:
            # Removing the -N suffix CellRanger adds
            labels.index = labels.index.str.replace(r"-\d+$", "", regex=True)
            common = adata.obs_names.intersection(labels.index)

        print(f"{donor}: {len(common)} barcodes matched between QC h5ad and Azimuth labels")
        adata = adata[common].copy()
        labels = labels.loc[common]

        # add Azimuth labels to obs
        adata.obs["predicted.celltype.l2"]       = labels["predicted.celltype.l2"].values
        adata.obs["predicted.celltype.l2.score"]  = labels["predicted.celltype.l2.score"].values
        adata.obs["predicted.celltype.l1"]       = labels["predicted.celltype.l1"].values
        adata.obs["donor"] = donor
        adata.obs["sex"]   = meta["sex"]
        adata.obs["age"]   = meta["age"]

        # filter to high confidence cells
        before = adata.n_obs
        adata = adata[adata.obs["predicted.celltype.l2.score"] >= min_score].copy()
        print(f"{donor}: {before} -> {adata.n_obs} cells after confidence filter (score >= {min_score})")

        # map to Memory_B / Plasma / Other
        adata.obs["cell_type"] = adata.obs["predicted.celltype.l2"].map(LABEL_MAP).fillna("Other")

        adatas.append(adata)

    if not adatas:
        raise ValueError("No donors loaded — check qc_dir and azimuth_dir paths")

    # combine all donors
    combined = adatas[0].concatenate(*adatas[1:], batch_key="donor_batch", index_unique="-")

    # cell type summary
    props = (combined.obs.groupby(["donor", "sex", "age", "cell_type"]).size()
             .reset_index(name="n_cells"))
    props["frac"] = props["n_cells"] / props.groupby("donor")["n_cells"].transform("sum")
    props.to_csv(f"{outdir}/celltype_proportions_by_donor.csv", index=False)
    print("\nCell type counts across cohort:")
    print(combined.obs["cell_type"].value_counts())

    # EDA plots
    sc.pl.umap(combined, color=["cell_type", "sex", "age", "donor"],
               save="_azimuth_celltype_sex_age_donor.png", show=False)
    sc.pl.umap(combined, color=["predicted.celltype.l2", "predicted.celltype.l2.score"],
               save="_azimuth_l2_labels.png", show=False)

    # write combined h5ad
    combined.write(f"{outdir}/combined_annotated.h5ad")
    print(f"\nCombined h5ad written: {outdir}/combined_annotated.h5ad")

    # write per donor x cell type h5ads for pySCENIC
    n_written = 0
    for donor in DONOR_META.keys():
        for ct in ["Memory_B", "Plasma"]:
            sub = combined[(combined.obs.donor == donor) &
                           (combined.obs.cell_type == ct)].copy()
            if sub.n_obs == 0:
                print(f"WARNING: {donor} {ct} has 0 cells after filtering")
                continue
            out_path = f"{outdir}/{donor}_{ct}.h5ad"
            sub.write(out_path)
            print(f"{donor} {ct}: {sub.n_obs} cells written -> {out_path}")
            n_written += 1

    print(f"\n{n_written}/16 donor x cell type h5ads written")
    print("Done.")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--qc_dir",      required=True)
    p.add_argument("--azimuth_dir", required=True)
    p.add_argument("--outdir",      required=True)
    p.add_argument("--min_score",   type=float, default=0.5)
    args = p.parse_args()
    main(args.qc_dir, args.azimuth_dir, args.outdir, args.min_score)
