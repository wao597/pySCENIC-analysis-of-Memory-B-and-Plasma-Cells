#!/usr/bin/env bash
#SBATCH --job-name=gprofiler
#SBATCH --cpus-per-task=60
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=gprofiler_out/logs/gprofiler_%j.out
#SBATCH --error=gprofiler_out/logs/gprofiler_%j.err
#
# 07_run_gprofiler.sh

#
# Usage:
#   mkdir -p gprofiler_out/logs
#   sbatch 07_run_gprofiler.sh


module purge
module load bluebear
module load bear-apps/2024a
module load scanpy/1.9.8-foss-2024a

set -euo pipefail

LASSO_DIR="/rds/projects/r/russdr-bb-data/wao597/lasso_out"
DE_CSV="/rds/projects/r/russdr-bb-data/wao597/pseudobulk_out/de_memoryB_vs_plasma.csv"
GRN_CSV="/rds/projects/r/russdr-bb-data/wao597/grn_analysis_out/most_common_pairs.csv"
OUTDIR="/rds/projects/r/russdr-bb-data/wao597/gprofiler_out"

mkdir -p "$OUTDIR"

# Confirming gprofiler is available
python -c "from gprofiler import GProfiler; print('gprofiler ok')"


python 07_gprofiler_enrichment.py \
    --lasso_dir "$LASSO_DIR" \
    --de_csv "$DE_CSV" \
    --grn_csv "$GRN_CSV" \
    --outdir "$OUTDIR" \
    --tfs BACH2 ATF4 XBP1 RORA\
    --n_de 200 \
    --top_terms 15


