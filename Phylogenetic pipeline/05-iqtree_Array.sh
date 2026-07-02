#!/bin/bash

#SBATCH -c 8
#SBATCH --mem=25G
#SBATCH -t 48:00:00
#SBATCH --array=0-4%5
#SBATCH --output=logs/iqtree_%A_%a.out
#SBATCH --error=logs/iqtree_%A_%a.err
#SBATCH --job-name=iqtree_array

# Activate IQ-TREE environment
source /path/to/miniconda3/bin/activate
source activate /path/to/apps/conda_envs/mamba/envs/iqtree

# Set paths
ALIGN_DIR="blast_region_LGT_check_Alignments"
TREE_DIR="blast_region_LGT_check_Trees"

mkdir -p "$TREE_DIR" logs

# Get list of alignment files
ALIGN_FILES=(${ALIGN_DIR}/*.fasta)
TOTAL=${#ALIGN_FILES[@]}

if [ $SLURM_ARRAY_TASK_ID -le $TOTAL ]; then
    alignment_file=${ALIGN_FILES[$((SLURM_ARRAY_TASK_ID-1))]}
    query_name=$(basename "$alignment_file" .fa)

    mkdir -p "${TREE_DIR}/${query_name}"

    echo "[$(date)] Building tree for: $query_name"

    iqtree -s "$alignment_file" \
           -m MFP \
           -B 1000 \
           -nt 8 \
           -pre "${TREE_DIR}/${query_name}/${query_name}"

    echo "[$(date)] Finished: $query_name"
else
    echo "No alignment file for task ID $SLURM_ARRAY_TASK_ID"
fi
