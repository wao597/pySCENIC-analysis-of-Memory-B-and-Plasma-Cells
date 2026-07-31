#!/usr/bin/env bash
#SBATCH --job-name=plasma_umap
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=plasma_umap_%j.out
#SBATCH --error=plasma_umap_%j.err

module purge
module load bluebear
module load bear-apps/2024a
module load scanpy/1.9.8-foss-2024a

python bach2_plasma_umap.py --h5ad /rds/projects/r/russdr-bb-data/wao597/pseudobulk_out/balanced_cells.h5ad --outdir /rds/projects/r/russdr-bb-data/wao597/bach2_coexpression_out
