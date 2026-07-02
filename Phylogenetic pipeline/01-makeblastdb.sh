#!/bin/bash

# Activate conda environment
source /path/to/miniconda3/bin/activate
source activate /path/to/miniconda3/envs/blast

# Make blast DBs of all chloroplasts

    makeblastdb -in path/to/chloroplast_assembly.fna \
    -dbtype nucl \
    -title [species] \
    -parse_seqids \
    -out path/to/BlastDB/[species]/[species]
