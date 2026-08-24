# Identification-and-phylogenetic-analysis-of-NUPTs
This pipeline includes the scripts used to identify and phylogenetically analyse NUPTs to infer their potential lateral acquisition.

The scripts are presented in numerical running order and the NUPT Identification Pipeline needs to be run prior to the Phylogenetic pipeline. 

### NUPT IDENTIFICATION PIPELINE ###

01-makeblastdb.sh
Generate BLAST databases for nuclear genomes.

02.1-self_blast.sh
BLAST nuclear genomes own chloroplast genome against itself.

02.2-Foreign_chloroplast_blast.sh
BLAST foreign chloroplasts against nuclear genome.

03-filter_blast_hits.sh
Filter the resulting self and foreign blast hits >=1kb & 95% Percentage identity

04-compare_chloroplast_blasts.py
Compare the filtered self and foreign blast hits to identify candidates where bitscore >100 for the foreign hit.

05-extract_foreign_NUPT_regions.sh
Extract candidate foreign hits as FASTA sequences.

### NUPT PHYLOGENETIC PIPELINE ###

01-makeblastdb.sh
Generate BLAST databases for chloroplast genomes.

02-blastn_loop.sh
BLAST candidate NUPTs against chloroplast BLAST databases.

03-filter_and_merge_1000bp_blasts.sh
Calls 03-extract_best_hsp.py to filter BLAST hits per chloroplast by length (>1kb), orientation (Highest cumulative bitscore & alignment length) and then best hit per chloroplast (Highest bitscore & alignment length).

04-align_blast_regions.sh
Align top hits per chloroplast using MAFFT.

05-iqtree_Array.sh
Construct maximum likelihood phylogenetic trees from the MAFFT MSA.

label_trees_v2.py
Label chloroplast accessions by there species names.
