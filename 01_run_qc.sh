#!/usr/bin/env bash

#SBATCH --job-name=qc

#SBATCH --array=0-7

#SBATCH --cpus-per-task=60

#SBATCH --mem=64G

#SBATCH --time=4:00:00


#SBATCH --output=qc_all_samples/logs/qc_%A_%a.out

#SBATCH --error=qc_all_samples/logs/qc_%A_%a.err

#

# 01_run_qc.sh

#

# SLURM array job: one task per donor, runs 01_qc.py for that donor.

#

# Usage:

#   mkdir -p qc_all_samples/logs

#   sbatch 01_run_qc.sh /path/to/cellranger_out /path/to/qc_all_samples



module purge
module load bluebear
module load bear-apps/2024a
module load scanpy/1.9.8-foss-2024a

set -euo pipefail



CELLRANGER_DIR=$1

OUTDIR=$2



# donor metadata 

DONORS=(BM1  BM2  BM3  BM4  BM5  BM6  BM7  BM8)

SEXES=(female male male male male female female female)

AGES=(52 50 39 29 29 26 36 32)



DONOR=${DONORS[$SLURM_ARRAY_TASK_ID]}

SEX=${SEXES[$SLURM_ARRAY_TASK_ID]}

AGE=${AGES[$SLURM_ARRAY_TASK_ID]}



mkdir -p qc_all_samples/logs

mkdir -p "$OUTDIR"



echo "=== Task $SLURM_ARRAY_TASK_ID: $DONOR (sex=$SEX, age=$AGE) ==="



python 01_qc.py --cellranger_dir "$CELLRANGER_DIR" --donor "$DONOR" --sex "$SEX" --age "$AGE" --outdir "$OUTDIR"


echo "=== $DONOR QC complete ==="
