#!/bin/bash

source /path/to/miniconda3/bin/activate
source activate /path/to/apps/conda_envs/mamba/envs/seqkit

module load Biopython

Sample_info="samples.txt"
QUERY_DIR="Blast_out_LGT_check"

FILTERED_HITS_DIR="${QUERY_DIR}/filtered_hits"
FASTA_OUT_DIR="${FILTERED_HITS_DIR}/blast_region_fasta"
BEST_HSP_DIR="${FILTERED_HITS_DIR}/best_hsp_fasta"

PYTHON_SCRIPT="extract_best_hsp.py"

mkdir -p "$FILTERED_HITS_DIR"
mkdir -p "$FASTA_OUT_DIR"
mkdir -p "$BEST_HSP_DIR"

for query_folder in ${QUERY_DIR}/*/; do
    query_name=$(basename "$query_folder")
    echo "Processing query: $query_name"

    mkdir -p "${FASTA_OUT_DIR}/${query_name}"
    mkdir -p "${BEST_HSP_DIR}/${query_name}"

    while read -r sample; do
        echo "  Processing sample: $sample"

        RAW_BLAST_FILE="${query_folder}/${sample}_${query_name}_blast.out"
        FILTERED_BLAST_FILE="${FILTERED_HITS_DIR}/${query_name}/${sample}_filtered_blast.out"
        FILTERED_FASTA="${FASTA_OUT_DIR}/${query_name}/${sample}_filtered_hsps.fasta"
        BEST_HSP_FASTA="${BEST_HSP_DIR}/${query_name}/${sample}_best_hsp.fasta"

        if [[ ! -f "$RAW_BLAST_FILE" ]]; then
            echo "  Warning: $RAW_BLAST_FILE not found. Skipping."
            continue
        fi

        echo "  Filtering hits (alignment length >= 1000)"
        mkdir -p "${FILTERED_HITS_DIR}/${query_name}"
        awk '$4 >= 1000' "$RAW_BLAST_FILE" > "$FILTERED_BLAST_FILE"

        if [[ ! -s "$FILTERED_BLAST_FILE" ]]; then
            echo "  No entries after filtering for $sample"
            continue
        fi

        echo "  Extracting HSP sequences"
        awk '{ printf(">%s_%s_%s\n%s\n", $2, $9, $10, $13) }' \
            "$FILTERED_BLAST_FILE" > "$FILTERED_FASTA"

        echo "  Running best HSP selection"
        python "$PYTHON_SCRIPT" "$FILTERED_BLAST_FILE" "$BEST_HSP_FASTA"

    done < "$Sample_info"

    # Optional: concatenate best HSP outputs per query
    cat "${BEST_HSP_DIR}/${query_name}"/*_best_hsp.fasta \
        > "${BEST_HSP_DIR}/${query_name}_best_hsp.fasta"

    # Optional: concatenate filtered HSPs per query
    cat "${FASTA_OUT_DIR}/${query_name}"/*_filtered_hsps.fasta \
        > "${FASTA_OUT_DIR}/${query_name}_filtered_hsps.fasta"

done

echo "All processing complete."
