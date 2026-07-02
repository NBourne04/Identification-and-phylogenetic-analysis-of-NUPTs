#!/usr/bin/env python3

#######################################################################
### This script is called from 03-filter_and_merge_1000bp_blasts.sh ###
#######################################################################

from collections import defaultdict
from Bio.Seq import Seq
import sys

if len(sys.argv) != 3:
    print("Usage: python extract_best_hsp.py input.tsv output.fasta")
    sys.exit(1)

blast_file = sys.argv[1]
output_fasta = sys.argv[2]

species_hits = defaultdict(list)

# -------------------------
# Read BLAST file
# -------------------------
with open(blast_file) as f:
    for line in f:
        fields = line.strip().split("\t")

        if len(fields) < 13:
            continue

        try:
            qstart = int(fields[6])
            qend   = int(fields[7])
            sstart = int(fields[8])
            send   = int(fields[9])
            bitscore = float(fields[11])
        except ValueError:
            continue

        seq = fields[12]
        sseqid = fields[1]

        q_start = min(qstart, qend)
        q_end   = max(qstart, qend)

        strand = "forward"
        if sstart > send:
            strand = "reverse"
            seq = str(Seq(seq).reverse_complement())

        species_hits[sseqid].append({
            "q_start": q_start,
            "q_end": q_end,
            "bitscore": bitscore,
            "length": q_end - q_start + 1,
            "strand": strand,
            "seq": seq
        })

# -------------------------
# Score orientation
# -------------------------
def score_hits(hits):
    total_bitscore = sum(h["bitscore"] for h in hits)
    total_length = sum(h["length"] for h in hits)
    return total_bitscore, total_length

# -------------------------
# Select best HSP within hits
# -------------------------
def select_best_hsp(hits):
    # Highest bitscore, tie-breaker = length
    return max(hits, key=lambda x: (x["bitscore"], x["length"]))

# -------------------------
# Process each accession
# -------------------------
with open(output_fasta, "w") as out:

    for species, hits in species_hits.items():

        forward_hits = [h for h in hits if h["strand"] == "forward"]
        reverse_hits = [h for h in hits if h["strand"] == "reverse"]

        if not forward_hits and not reverse_hits:
            continue

        f_score, f_len = score_hits(forward_hits) if forward_hits else (0, 0)
        r_score, r_len = score_hits(reverse_hits) if reverse_hits else (0, 0)

        # -------------------------
        # Choose best orientation
        # -------------------------
        if f_score > r_score:
            best_hits = forward_hits
            strand = "forward"

        elif r_score > f_score:
            best_hits = reverse_hits
            strand = "reverse"

        else:
            if f_len > r_len:
                best_hits = forward_hits
                strand = "forward"
            else:
                best_hits = reverse_hits
                strand = "reverse"

        if not best_hits:
            continue

        # -------------------------
        # Select best HSP within orientation
        # -------------------------
        best_hsp = select_best_hsp(best_hits)

        # -------------------------
        # Output
        # -------------------------
        header = f">{species}|{strand}|{best_hsp['q_start']}-{best_hsp['q_end']}|bitscore={best_hsp['bitscore']}"
        out.write(f"{header}\n{best_hsp['seq']}\n")

print("Done.")
