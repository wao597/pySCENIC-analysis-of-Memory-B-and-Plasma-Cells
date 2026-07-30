#!/usr/bin/env python
"""
03_h5ad_to_loom.py
Converts h5ad to loom for pySCENIC, filtering out genes without proper
HGNC symbols (Ensembl IDs, LINC genes starting with ENSG) so that
cisTarget can map them to the ranking database.

Usage:
    python 03_h5ad_to_loom.py /path/to/input.h5ad /path/to/output.loom
"""
import sys
import numpy as np
import h5py
import loompy
import scipy.sparse as sp

h5ad_path = sys.argv[1]
loom_path = sys.argv[2]

with h5py.File(h5ad_path, "r") as f:
    # read gene names
    var_index = f["var"]["_index"][:]
    gene_names = np.array([g.decode() if isinstance(g, bytes) else g for g in var_index])

    # read cell barcodes
    obs_index = f["obs"]["_index"][:]
    cell_ids = np.array([c.decode() if isinstance(c, bytes) else c for c in obs_index])

    # read raw counts layer
    counts_group = f["layers"]["counts"]
    if "data" in counts_group:
        data = counts_group["data"][:]
        indices = counts_group["indices"][:]
        indptr = counts_group["indptr"][:]
        shape = (len(cell_ids), len(gene_names))
        X = sp.csr_matrix((data, indices, indptr), shape=shape)
        X_dense = X.toarray()
    else:
        X_dense = counts_group[:]

# filter to HGNC symbols only 
hgnc_mask = np.array([not g.startswith("ENSG") for g in gene_names])
gene_names_filtered = gene_names[hgnc_mask]
X_filtered = X_dense[:, hgnc_mask]

print(f"Genes before filtering: {len(gene_names)}")
print(f"Genes after filtering (HGNC only): {len(gene_names_filtered)}")
print(f"Removed: {len(gene_names) - len(gene_names_filtered)} Ensembl-ID genes")

row_attrs = {"Gene": gene_names_filtered}
col_attrs = {"CellID": cell_ids}

loompy.create(loom_path, X_filtered.T, row_attrs, col_attrs)
print(f"Written: {loom_path} ({len(cell_ids)} cells x {len(gene_names_filtered)} genes)")
