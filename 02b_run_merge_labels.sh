#!/usr/bin/env bash
#SBATCH --job-name=merge_labels
#SBATCH --cpus-per-task=70
#SBATCH --mem=256G
#SBATCH --time=4:00:00
#SBATCH --output=azimuth_out/logs/merge_%j.out
#SBATCH --error=azimuth_out/logs/merge_%j.err
#
# 02b_run_merge_labels.sh
#
# This script merges Azimuth labels into QC h5ads, subsets to Memory B / Plasma, writes combined h5ad and per-donor/cell-type h5ads.

#
# Usage:
#   sbatch 02b_run_merge_labels.sh

module purge
module load bluebear
module load bear-apps/2024a
module load scanpy/1.9.8-foss-2024a

set -euo pipefail

QC_DIR="/rds/projects/r/russdr-bb-data/wao597/qc_all_samples"
AZIMUTH_DIR="/rds/projects/r/russdr-bb-data/wao597/azimuth_out"
OUTDIR="/rds/projects/r/russdr-bb-data/wao597/celltype_out"

mkdir -p "$OUTDIR"



python 02b_merge_azimuth_labels.py --qc_dir "$QC_DIR" --azimuth_dir "$AZIMUTH_DIR" --outdir "$OUTDIR" --min_score 0.5


