#!/usr/bin/env bash
#SBATCH --job-name=pyscenic
#SBATCH --array=0-15
#SBATCH --cpus-per-task=70
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=pyscenic_out/logs/pyscenic_%A_%a.out
#SBATCH --error=pyscenic_out/logs/pyscenic_%A_%a.err
#
# 03_run_pyscenic.sh
#
# This is an array array job: 16 tasks total (8 donors x 2 cell types).
# Array index mapping:
#   0-7  = Memory_B (BM1-BM8)
#   8-15 = Plasma   (BM1-BM8)
#

#
# Usage:
#   mkdir -p pyscenic_out/logs
#   sbatch 03_run_pyscenic.sh

module purge
module load bluebear
module load bear-apps/2024a
module load pySCENIC/0.12.1-20250109-foss-2024a
export PYTHONNOUSERSITE=1

set -euo pipefail

#hardcoded paths
CELLTYPE_DIR="/rds/projects/r/russdr-bb-data/wao597/celltype_out"
OUTDIR="/rds/projects/r/russdr-bb-data/wao597/pyscenic_out"
TF_LIST="/rds/projects/r/russdr-bb-data/wao597/data/hs_hgnc_curated_tfs.txt"
RANKING_DB1="/rds/projects/r/russdr-bb-data/wao597/data/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
RANKING_DB2="/rds/projects/r/russdr-bb-data/wao597/data/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
MOTIF_TBL="/rds/projects/r/russdr-bb-data/wao597/data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

#array index - donor & cell type mapping
DONORS=(BM1 BM2 BM3 BM4 BM5 BM6 BM7 BM8 BM1 BM2 BM3 BM4 BM5 BM6 BM7 BM8)
CELLTYPES=(Memory_B Memory_B Memory_B Memory_B Memory_B Memory_B Memory_B Memory_B Plasma Plasma Plasma Plasma Plasma Plasma Plasma Plasma)

DONOR=${DONORS[$SLURM_ARRAY_TASK_ID]}
CELLTYPE=${CELLTYPES[$SLURM_ARRAY_TASK_ID]}
NAME="${DONOR}_${CELLTYPE}"

H5AD="${CELLTYPE_DIR}/${NAME}.h5ad"
RUNDIR="${OUTDIR}/${NAME}"

mkdir -p "$RUNDIR"
mkdir -p pyscenic_out/logs

echo "=== Task $SLURM_ARRAY_TASK_ID: $NAME ==="

# Double checking input exists
if [[ ! -f "$H5AD" ]]; then
    echo "ERROR: input not found: $H5AD" >&2
    exit 1
fi

# 1. Convert h5ad to loom
echo "=== $NAME: converting to loom ==="
python 03_h5ad_to_loom.py "$H5AD" "$RUNDIR/expr.loom"

#2. GRN inference (GRNBoost2)
pyscenic grn "$RUNDIR/expr.loom" "$TF_LIST" -o "$RUNDIR/adjacencies.tsv" --num_workers ${SLURM_CPUS_PER_TASK:-16}

#3: cisTarget pruning 
pyscenic ctx "$RUNDIR/adjacencies.tsv" "$RANKING_DB1" "$RANKING_DB2" --annotations_fname "$MOTIF_TBL" --expression_mtx_fname "$RUNDIR/expr.loom" --output "$RUNDIR/regulons.csv" --num_workers ${SLURM_CPUS_PER_TASK:-16}

#4: AUCell scoring 
pyscenic aucell "$RUNDIR/expr.loom" "$RUNDIR/regulons.csv" -o "$RUNDIR/auc_mtx.csv" --num_workers ${SLURM_CPUS_PER_TASK:-16}


