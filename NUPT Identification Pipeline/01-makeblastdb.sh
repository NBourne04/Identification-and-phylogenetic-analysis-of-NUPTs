#!/bin/bash

# Activate conda environment
source /path/to/miniconda3/bin/activate
source activate /path/to/miniconda3/envs/blast

    makeblastdb -in path/to/Genome_assembly.fna \
    -dbtype nucl \
    -title [species] \
    -parse_seqids \
    -out path/to/BlastDB/[species]/[species]
