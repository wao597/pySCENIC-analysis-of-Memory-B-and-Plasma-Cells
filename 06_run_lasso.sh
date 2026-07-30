#!/usr/bin/env bash
#SBATCH --job-name=lasso_tf
#SBATCH --cpus-per-task=60
#SBATCH --mem=128G
#SBATCH --time=6:00:00
#SBATCH --output=lasso_out/%x_%j.out
#SBATCH --error=lasso_out/%x_%j.err
#
# 06_run_lasso.sh
#
# Trains a Lasso model per cell type for a TF of interest.
#
# To run with the default TF (PRDM1):
#   mkdir -p lasso_out
#   sbatch 06_run_lasso.sh
#
# To swap to a different TF, mention it as an argument:
#   sbatch 06_run_lasso.sh BACH2
#   sbatch 06_run_lasso.sh ATF4
#   sbatch 06_run_lasso.sh RORA
#


module purge
module load bluebear
module load bear-apps/2024a
module load scanpy/1.9.8-foss-2024a

set -euo pipefail

# PRDM1 is the default TF bit specify if investigating another
TF=${1:-PRDM1}

BALANCED_H5AD="/rds/projects/r/russdr-bb-data/wao597/pseudobulk_out/balanced_cells.h5ad"
OUTDIR="/rds/projects/r/russdr-bb-data/wao597/lasso_out/${TF}"

mkdir -p "$OUTDIR"
mkdir -p lasso_out

echo "=== Lasso model for TF: $TF ==="
echo "=== Output dir: $OUTDIR ==="

python 06_lasso_tf_model.py --balanced_h5ad "$BALANCED_H5AD" --tf "$TF" --outdir "$OUTDIR"

echo "=== Lasso complete for $TF ==="
