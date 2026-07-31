#!/usr/bin/env bash
#SBATCH --job-name=celltype
#SBATCH --cpus-per-task=60
#SBATCH --mem=256G
#SBATCH --time=4:00:00
#SBATCH --output=celltype_out/logs/celltype_%j.out
#SBATCH --error=celltype_out/logs/celltype_%j.err
#
# 02_run_celltype.sh

# Usage:
#   mkdir -p celltype_out/logs
#   sbatch 02_run_celltype.sh /path/to/qc_all_samples /path/to/celltype_out

module purge
module load bluebear
module load bear-apps/2024a
module load scanpy/1.9.8-foss-2024a

set -euo pipefail

QC_DIR=$1
OUTDIR=$2

mkdir -p "$OUTDIR"

python 02_celltype_annotation.py --qc_dir "$QC_DIR" --outdir "$OUTDIR"

