#!/usr/bin/env bash
#SBATCH --job-name=overlap_analysis
#SBATCH --cpus-per-task=60
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=overlap_out/logs/overlap_%j.out
#SBATCH --error=overlap_out/logs/overlap_%j.err
#
# 10_run_overlap.sh
# This script overlap analysis connecting Memory B and Plasma cell findings.
#
# Usage:
#   mkdir -p overlap_out/logs
#   sbatch 10_run_overlap.sh

module purge
module load bluebear
module load bear-apps/2024a
module load scanpy/1.9.8-foss-2024a

set -euo pipefail

PYSCENIC_DIR="/rds/projects/r/russdr-bb-data/wao597/pyscenic_out"
LASSO_DIR="/rds/projects/r/russdr-bb-data/wao597/lasso_out"
DE_CSV="/rds/projects/r/russdr-bb-data/wao597/pseudobulk_out/de_memoryB_vs_plasma.csv"
PAIRS_CSV="/rds/projects/r/russdr-bb-data/wao597/grn_analysis_out/most_common_pairs.csv"
OUTDIR="/rds/projects/r/russdr-bb-data/wao597/overlap_out"

mkdir -p "$OUTDIR"


python 10_overlap_analysis.py \
    --pyscenic_dir "$PYSCENIC_DIR" \
    --lasso_dir    "$LASSO_DIR" \
    --de_csv       "$DE_CSV" \
    --pairs_csv    "$PAIRS_CSV" \
    --outdir       "$OUTDIR" \
    --tfs BACH2 ATF4 XBP1


