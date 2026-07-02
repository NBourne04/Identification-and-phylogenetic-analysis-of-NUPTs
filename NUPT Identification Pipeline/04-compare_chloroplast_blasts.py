#!/usr/bin/env python3
"""
04-compare_chloroplast_blasts.py

Purpose:
Compare whether any foreign grass chloroplast genome has a BETTER hit
than Alloteropsis semialata chloroplast to the SAME nuclear region.

Usage: 04-compare_chloroplast_blasts.py self_blast foreign_blast

FORMAT REQUIRED (outfmt 6):
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

"Better" hit ranking:
1) Higher bit score
2) Higher percent identity

FOREIGN is only called if:
bitscore difference >=100

OUTPUT:
One line per self chloroplast hit, showing best competing foreign hit (if any)
"""

import argparse
from collections import defaultdict


def parse_blast(file):
    hits = defaultdict(list)

    with open(file) as f:
        for line in f:
            if not line.strip():
                continue
            p = line.strip().split("\t")

            hit = {
                "qseqid": p[0],
                "sseqid": p[1],
                "pident": float(p[2]),
                "length": int(p[3]),
                "qstart": int(p[6]),
                "qend": int(p[7]),
                "sstart": int(p[8]),
                "send": int(p[9]),
                "evalue": float(p[10]),
                "bitscore": float(p[11]),
            }

            hits[p[1]].append(hit)

    return hits


def overlap(a1, a2, b1, b2):
    start = max(min(a1, a2), min(b1, b2))
    end = min(max(a1, a2), max(b1, b2))
    return max(0, end - start + 1)


def species_from_query(qid):
    if qid.endswith("_chloroplast"):
        return qid.rsplit("_chloroplast", 1)[0]
    else:
        return qid


def better_hit(h1, h2):
    """Return True if h1 is better than h2"""
    if h1["bitscore"] != h2["bitscore"]:
        return h1["bitscore"] > h2["bitscore"]
    return h1["pident"] > h2["pident"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("self_blast")
    ap.add_argument("foreign_blast")
    ap.add_argument("--min_overlap", type=int, default=50)
    ap.add_argument("-o", default="NUPT_competition.tsv")
    args = ap.parse_args()

    self_hits = parse_blast(args.self_blast)
    foreign_hits = parse_blast(args.foreign_blast)

    out = open(args.o, "w")

    header = [
        "Nuclear_seq",
        "Nuc_start", "Nuc_end",
        "Self_bitscore", "Self_pident", "Self_length",
        "Best_foreign_species",
        "Foreign_bitscore", "Foreign_pident", "Foreign_length",
        "Bitscore_diff",
        "Winner"
    ]
    out.write("\t".join(header) + "\n")

    total = 0
    foreign_wins = 0

    for nuc in self_hits:
        if nuc not in foreign_hits:
            continue

        for self_hit in self_hits[nuc]:

            best_foreign = None

            for f in foreign_hits[nuc]:
                ov = overlap(self_hit["sstart"], self_hit["send"],
                             f["sstart"], f["send"])
                if ov < args.min_overlap:
                    continue

                if best_foreign is None or better_hit(f, best_foreign):
                    best_foreign = f

            if best_foreign is None:
                continue

            total += 1

            delta = best_foreign["bitscore"] - self_hit["bitscore"]

            if delta >= 100:
                winner = "FOREIGN"
                foreign_wins += 1
            else:
                winner = "SELF"

            row = [
                nuc,
                str(self_hit["sstart"]), str(self_hit["send"]),
                f"{self_hit['bitscore']}", f"{self_hit['pident']}", str(self_hit["length"]),
                species_from_query(best_foreign["qseqid"]),
                f"{best_foreign['bitscore']}", f"{best_foreign['pident']}", str(best_foreign["length"]),
                f"{delta}",
                winner
            ]

            out.write("\t".join(row) + "\n")

    out.close()

    print(f"Self hits evaluated: {total}")
    print(f"Foreign better than self (Δbitscore ≥ 20): {foreign_wins}")


if __name__ == '__main__':
    main()
