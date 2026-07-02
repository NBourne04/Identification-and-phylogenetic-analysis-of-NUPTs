#!/bin/bash

# Activate conda environment
source /path/to/miniconda3/bin/activate
source activate /path/to/miniconda3/envs/blast

# List of species
species_list=("AUS1_2019" "RSA5-3" "TAN1-04B" "ZAM1505-10")

for species in "${species_list[@]}"; do
    echo "Running BLAST for $species..."

    blastn -query Chloroplasts/NC_027824.1.fasta \
           -db "BlastDB/${species}/${species}" \
           -out "Blast_outs/ASEM/ASEM_chloroplast_vs_${species}.out" \
           -evalue 1e-10 \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" \
           -num_threads 1

done
