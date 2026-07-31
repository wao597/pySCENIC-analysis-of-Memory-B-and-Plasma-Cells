#!/usr/bin/env bash
#SBATCH --job-name=pseudobulk_de
#SBATCH --cpus-per-task=60
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=pseudobulk_out/logs/pseudobulk_%j.out
#SBATCH --error=pseudobulk_out/logs/pseudobulk_%j.err
#
# 05_run_pseudobulk_de.sh
#
# This is a single job, balanced cell sampling per donor per cell type, pseudobulk aggregation, and DESeq2 DE between Memory B and Plasma cells.
#
# Usage:
#   mkdir -p pseudobulk_out/logs
#   sbatch 05_run_pseudobulk_de.sh

module purge
module load bluebear
module load bear-apps/2024a
module load scanpy/1.9.8-foss-2024a

set -euo pipefail

COMBINED_H5AD="/rds/projects/r/russdr-bb-data/wao597/celltype_out/combined_annotated.h5ad"
OUTDIR="/rds/projects/r/russdr-bb-data/wao597/pseudobulk_out"

mkdir -p "$OUTDIR"

python -c "from pydeseq2.dds import DeseqDataSet; print('PyDESeq2 ok')"

python 05_pseudobulk_de.py --combined_h5ad "$COMBINED_H5AD" --outdir "$OUTDIR" --seed 0


