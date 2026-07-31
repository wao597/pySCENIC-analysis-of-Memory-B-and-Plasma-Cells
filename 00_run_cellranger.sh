#!/usr/bin/env bash
#SBATCH --job-name=cellranger
#SBATCH --array=0-7
#SBATCH --cpus-per-task=70
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/cellranger_%A_%a.out
#SBATCH --error=logs/cellranger_%A_%a.err

# This script finds fastq "sample" names per donor inside --fastqs dir and runs a cellranger count job per donor combining all of that donor's runs.


# Usage:
#   ./00_run_cellranger.sh /path/to/all_samples /path/to/refdata-gex-GRCh38-2024-A /path/to/outdir


module purge
module load bluebear
module load CellRanger/7.0.0



set -euo pipefail



FASTQ_DIR=$1

TRANSCRIPTOME=$2

OUTDIR=$3

DONORS=(BM1 BM2 BM3 BM4 BM5 BM6 BM7 BM8)



DONOR=${DONORS[$SLURM_ARRAY_TASK_ID]}

DONOR_FASTQ_DIR="$FASTQ_DIR/Manton${DONOR}"



mkdir -p "$OUTDIR"

cd "$OUTDIR"



if [[ ! -d "$DONOR_FASTQ_DIR" ]]; then

    echo "ERROR: directory not found: $DONOR_FASTQ_DIR" >&2

    exit 1

fi



SAMPLES=$(ls "$DONOR_FASTQ_DIR" | grep -E "^Manton${DONOR}_" | sed -E 's/_S[0-9]+_L[0-9]+_R[12]_001\.fastq\.gz$//' | sort -u | paste -sd, -)



if [[ -z "$SAMPLES" ]]; then

    echo "ERROR: no fastqs found for $DONOR in $DONOR_FASTQ_DIR" >&2

    exit 1

fi



echo "=== Task $SLURM_ARRAY_TASK_ID: $DONOR -> samples: $SAMPLES ==="


cellranger count --id="Manton${DONOR}" --transcriptome="$TRANSCRIPTOME" --fastqs="$DONOR_FASTQ_DIR" --sample="$SAMPLES" --localcores=${SLURM_CPUS_PER_TASK:-16} --localmem=120

