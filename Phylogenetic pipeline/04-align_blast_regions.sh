#!/bin/bash

# Activate MAFFT environment
source /path/to/miniconda3/bin/activate
source activate /path/to/apps/conda_envs/mamba/envs/mafftv7.453

# Set up paths
FASTA_OUT_DIR="Blast_out_LGT_check/filtered_hits/best_hsp_fasta"
ALIGN_DIR="blast_region_LGT_check_Alignments"
Queries_info="queries/queries.txt"
QUERY_DIR="queries"

mkdir -p "$ALIGN_DIR"

# Loop through each query fasta filename (including .fasta) from queries.txt
while read -r query_fasta_file; do
    echo "Processing query file: $query_fasta_file"

    query_path="${QUERY_DIR}/${query_fasta_file}"
    query_name="${query_fasta_file%.fasta}"
    fragments_fasta="${FASTA_OUT_DIR}/${query_name}_best_hsp.fasta"
    output_alignment="${ALIGN_DIR}/aligned_${query_name}_best_hsp.fasta"

    # Check both files exist
    if [[ -f "$query_path" && -f "$fragments_fasta" ]]; then
        echo "  Aligning fragments to query..."

        mafft --auto \
              --adjustdirectionaccurately \
              --thread 8 \
              --maxiterate 1000 \
              --addfragments "$fragments_fasta" "$query_path" > "$output_alignment"
    else
        echo "  Warning: Missing query or fragments for $query_name"
    fi

done < "$Queries_info"

echo "All alignments complete."
