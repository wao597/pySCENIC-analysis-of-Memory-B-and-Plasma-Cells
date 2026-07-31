#!/usr/bin/env bash
#SBATCH --job-name=grn_analysis
#SBATCH --cpus-per-task=60
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=grn_analysis_out/logs/grn_%j.out
#SBATCH --error=grn_analysis_out/logs/grn_%j.err
#
# 04_run_grn_analysis.sh

# This job aggregates GRNBoost2 importance scores across all 16
# pySCENIC runs and generates violin+box plots of most recurrent pairs.
#
# Usage:
#   mkdir -p grn_analysis_out/logs
#   sbatch 04_run_grn_analysis.sh

module purge
module load bluebear
module load bear-apps/2024a
module load pySCENIC/0.12.1-20250109-foss-2024a
module load scanpy/1.9.8-foss-2024a

set -euo pipefail

PYSCENIC_DIR="/rds/projects/r/russdr-bb-data/wao597/pyscenic_out"
OUTDIR="/rds/projects/r/russdr-bb-data/wao597/grn_analysis_out"

mkdir -p "$OUTDIR"



python 04_aggregate_grn_scores.py --pyscenic_dir "$PYSCENIC_DIR" --outdir "$OUTDIR" --top_n 20


