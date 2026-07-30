#!/usr/bin/env bash
#SBATCH --job-name=tf_network
#SBATCH --cpus-per-task=60
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=network_out/logs/network_%j.out
#SBATCH --error=network_out/logs/network_%j.err
#
# 08_run_network_focused.sh
#
# This script builds BACH2 (Memory B) and ATF4 (Plasma) regulon network visual

#
# Usage:
#   mkdir -p network_out/logs
#   sbatch 08_run_network_focused.sh

module purge
module load bluebear
module load bear-apps/2024a
module load scanpy/1.9.8-foss-2024a

python -c "import networkx" 2>/dev/null || pip install --user networkx

set -euo pipefail

PAIRS_CSV="/rds/projects/r/russdr-bb-data/wao597/grn_analysis_out/most_common_pairs.csv"
ADJ_DIR="/rds/projects/r/russdr-bb-data/wao597/pyscenic_out"
OUTDIR="/rds/projects/r/russdr-bb-data/wao597/network_out"

mkdir -p "$OUTDIR"

echo "=== Building BACH2 and ATF4 regulon networks ==="

python 08_tf_network_focused.py \
    --pairs_csv "$PAIRS_CSV" \
    --adj_dir   "$ADJ_DIR" \
    --outdir    "$OUTDIR" \
    --top_targets 20 \
    --min_donors 6


