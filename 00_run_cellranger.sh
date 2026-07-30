#!/usr/bin/env bash
# 00_run_cellranger.sh
#
# This script finds fastq "sample" names per donor inside --fastqs dir
# (e.g. MantonBM1_HiSeq_1, MantonBM1_HiSeq_5, ...) and runs ONE cellranger
# count job per donor combining all of that donor's runs.
#
# Usage:
#   ./00_run_cellranger.sh /path/to/all_samples /path/to/refdata-gex-GRCh38-2024-A /path/to/outdir

set -euo pipefail

FASTQ_DIR=$1
TRANSCRIPTOME=$2
OUTDIR=$3
DONORS=(BM1 BM2 BM3 BM4 BM5 BM6 BM7 BM8)

mkdir -p "$OUTDIR"
cd "$OUTDIR"

for DONOR in "${DONORS[@]}"; do
    # Find every fastq "sample" stem belonging to this donor.
    # cellranger sample names are the bit before _S*_L*_R*_001.fastq.gz,
    # e.g. MantonBM1_HiSeq_1_S1_L001_R1_001.fastq.gz -> MantonBM1_HiSeq_1
    SAMPLES=$(ls "$FASTQ_DIR" \
        | grep -E "^Manton${DONOR}_" \
        | sed -E 's/_S[0-9]+_L[0-9]+_R[12]_001\.fastq\.gz$//' \
        | sort -u \
        | paste -sd, -)

    if [[ -z "$SAMPLES" ]]; then
        echo "WARNING: no fastqs found for $DONOR, skipping" >&2
        continue
    fi

    echo "=== $DONOR -> samples: $SAMPLES ==="

    cellranger count \
        --id="Manton${DONOR}" \
        --transcriptome="$TRANSCRIPTOME" \
        --fastqs="$FASTQ_DIR" \
        --sample="$SAMPLES" \
        --create-bam=true \
        --localcores=16 \
        --localmem=128
done
