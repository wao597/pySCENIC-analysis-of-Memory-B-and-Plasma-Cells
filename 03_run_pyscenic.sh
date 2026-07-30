#!/usr/bin/env bash
# 03_run_pyscenic.sh
#
# Runs the 3-stage pySCENIC workflow (GRNBoost2 -> cisTarget -> AUCell)
# separately for every donor x cell_type combination (16 runs for 8 donors x
# {Memory_B, Plasma}), exactly as requested ("run pySCENIC on each individually").
#
# Requires: pyscenic installed, cisTarget feather rankings + motif annotation
# tsv downloaded (https://resources.aertslab.org/cistarget/), and a TF list
# (e.g. allTFs_hg38.txt from the same resource).
#
# Usage:
#   ./03_run_pyscenic.sh ./celltype_out ./pyscenic_out \
#       /path/to/hg38_TFs.txt \
#       /path/to/hg38__refseq-r80__*.feather \
#       /path/to/motifs-v10-nr.hgnc-m0.001-o0.0.tbl

set -euo pipefail
INDIR=$1
OUTDIR=$2
TF_LIST=$3
RANKING_DB=$4   # can be a glob, quote it
MOTIF_TBL=$5

mkdir -p "$OUTDIR"

for H5AD in "$INDIR"/*_Memory_B.h5ad "$INDIR"/*_Plasma.h5ad; do
    NAME=$(basename "$H5AD" .h5ad)     # e.g. BM1_Memory_B
    RUNDIR="$OUTDIR/$NAME"
    mkdir -p "$RUNDIR"

    # pyscenic wants a loom file with raw counts
    python - "$H5AD" "$RUNDIR/expr.loom" <<'PYEOF'
import sys, scanpy as sc, loompy, numpy as np
adata = sc.read_h5ad(sys.argv[1])
X = adata.layers["counts"]
row_attrs = {"Gene": np.array(adata.var_names)}
col_attrs = {"CellID": np.array(adata.obs_names)}
loompy.create(sys.argv[2], X.T.toarray() if hasattr(X, "toarray") else X.T,
              row_attrs, col_attrs)
PYEOF

    echo "=== $NAME: GRN inference ==="
    pyscenic grn "$RUNDIR/expr.loom" "$TF_LIST" \
        -o "$RUNDIR/adjacencies.tsv" \
        --num_workers 8

    echo "=== $NAME: cisTarget pruning ==="
    pyscenic ctx "$RUNDIR/adjacencies.tsv" $RANKING_DB \
        --annotations_fname "$MOTIF_TBL" \
        --expression_mtx_fname "$RUNDIR/expr.loom" \
        --output "$RUNDIR/regulons.csv" \
        --num_workers 8

    echo "=== $NAME: AUCell scoring ==="
    pyscenic aucell "$RUNDIR/expr.loom" "$RUNDIR/regulons.csv" \
        -o "$RUNDIR/auc_mtx.csv" \
        --num_workers 8
done
