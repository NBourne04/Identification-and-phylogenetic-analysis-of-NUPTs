#!/bin/bash

# Activate conda environment
source /path/to/miniconda3/bin/activate
source activate /path/to/miniconda3/envs/blast

mkdir -p Blast_out_LGT_check

# Sample and query file lists

# text file listing the chloroplast accesions
Sample_info=samples.txt

# text file listing the NUPT candidate fasta files
Queries_info=queries/queries.txt

# Loop through all queries
while read -r query; do
    query_path="queries/${query}"
    query_name=$(basename "$query" .fasta)  # Adjust extension if needed

    # Create a subdirectory in Blast_out for this query
    mkdir -p "Blast_out_LGT_check/${query_name}"

    # Loop through all samples
    while read -r sample; do
        echo "Running BLAST for query: ${query} against sample: ${sample}"

        blastn -query "$query_path" \
               -db "../BlastDB/${sample}/${sample}" \
               -out "Blast_out_LGT_check/${query_name}/${sample}_${query_name}_blast.out" \
               -evalue 1e-10 \
               -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" \
               -num_threads 1
    done < "$Sample_info"
done < "$Queries_info"
