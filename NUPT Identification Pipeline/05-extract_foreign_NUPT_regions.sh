#!/bin/bash

module load BEDTools/2.31.0-GCC-12.3.0

set -euo pipefail

### AUS1_2019 RSA5-3 TAN1-04B ZAM1505-10 ###

species=ZAM1505-10

TSV=NUPT_blast_comparisons/${species}_NUPT_competition.tsv
GENOME=../Genomes/${species}/GCA_036785565.1_ZAM1505-10_genomic.fna
OUT=NUPT_blast_comparisons/${species}_NUPT_LGT_candidates.fasta

TMP_BED=$(mktemp)
TMP_SORTED=$(mktemp)
TMP_MERGED=$(mktemp)

echo "Filtering FOREIGN hits and converting to BED..."

awk 'NR>1 && $NF=="FOREIGN" {

    # Convert BLAST (1-based) BED (0-based)
    start = ($2 < $3) ? $2-1 : $3-1
    end   = ($2 < $3) ? $3   : $2
    strand = ($2 < $3) ? "+" : "-"
 
    print $1"\t"start"\t"end"\t.\t0\t"strand
}' "$TSV" | sed -E 's/gb\|([^|]+)\|/\1/' > "$TMP_BED"

echo "DEBUG cleaned BED:"
head -n 5 "$TMP_BED"

echo "Unique chromosomes:"
cut -f1 "$TMP_BED" | sort | uniq | head

echo "Sorting BED..."

bedtools sort -i "$TMP_BED" > "$TMP_SORTED"

echo "Merging overlapping regions..."

# Merge overlapping intervals (strand-aware merge is NOT default,
# so we temporarily ignore strand during merging for simplicity)
bedtools merge -i "$TMP_SORTED" > "$TMP_MERGED"

echo "Extracting merged sequences..."

bedtools getfasta \
    -fi "$GENOME" \
    -bed "$TMP_MERGED" \
    -fo "$OUT"

rm "$TMP_BED" "$TMP_SORTED" "$TMP_MERGED"

echo "Done. Output written to $OUT"
