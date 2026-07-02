#!/bin/bash

# Filter blast outputs for 95% percentage identity and 1000bp

# List of species
species_list=("AUS1_2019" "RSA5-3" "TAN1-04B" "ZAM1505-10")

for species in "${species_list[@]}"; do
    echo "Processing $species..."

    awk '$3 >= 95 && $4 >= 1000 { for (i=1; i<NF; i++) printf "%s\t", $i; printf "\n" }' \
        Blast_outs/ASEM/ASEM_chloroplast_vs_${species}.out \
        > Blast_outs/ASEM/ASEM_chloroplast_vs_${species}_filtered.out

    awk '$3 >= 95 && $4 >= 1000 { for (i=1; i<NF; i++) printf "%s\t", $i; printf "\n" }' \
        Blast_outs/outgroups/All_chloroplasts_no_ASEM_vs_${species}.out \
        > Blast_outs/outgroups/All_chloroplasts_no_ASEM_vs_${species}_filtered.out

done
