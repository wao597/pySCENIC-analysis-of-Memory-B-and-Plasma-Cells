#!/usr/bin/env bash
#SBATCH --job-name=azimuth
#SBATCH --array=0-7
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=128G
#SBATCH --time=01:00:00
#SBATCH --output=azimuth_out/logs/azimuth_%A_%a.out
#SBATCH --error=azimuth_out/logs/azimuth_%A_%a.err
#
# 02a_run_azimuth.sh
#
# Array job: one task per donor (8 donors -> --array=0-7)

# To use:
#   mkdir -p azimuth_out/logs
#   sbatch 02a_run_azimuth.sh

module purge
module load bluebear
module load bear-apps/2024a
module load R/4.5.0-gfbf-2024a
module load R-bundle-Bioconductor/3.21-foss-2024a-R-4.5.0

set -euo pipefail

# Defining paths
CELLRANGER_DIR="/rds/projects/r/russdr-bb-data/wao597/cellranger_out"
OUTDIR="/rds/projects/r/russdr-bb-data/wao597/azimuth_out"

# Setting up donor array
DONORS=(BM1 BM2 BM3 BM4 BM5 BM6 BM7 BM8)
DONOR=${DONORS[$SLURM_ARRAY_TASK_ID]}

mkdir -p "$OUTDIR"
mkdir -p azimuth_out/logs

echo "=== Task $SLURM_ARRAY_TASK_ID: $DONOR ==="

Rscript 02a_prepare_and_azimuth.R "$DONOR" "$CELLRANGER_DIR" "$OUTDIR"

echo "=== $DONOR Azimuth complete ==="
