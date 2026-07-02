#!/usr/bin/env python3

from Bio import Phylo
import re
import os
from pathlib import Path

# ---------- INPUT ----------

## species.txt needs to be the following tsv format: 
#      Accession    Species    Family    Subfamily
# e.g. NC_034680.1	Cynodon dactylon	Chloridoideae	Cynodonteae 

species_file = "species.txt"
output_dir = Path("labelled_trees")

output_dir.mkdir(exist_ok=True)

# ---------- BUILD ACCESSION MAP ----------
acc_to_species = {}

with open(species_file) as f:
    header = next(f)  # skip header line

    for line in f:
        parts = line.strip().split("\t")

        if len(parts) < 2:
            continue

        acc = parts[0]

        # Species → Genus_species
        species_name = parts[1].replace(" ", "_")

        # Columns (handle missing safely)
        subfamily = parts[2] if len(parts) > 2 and parts[2] else ""
        tribe = parts[3] if len(parts) > 3 and parts[3] else ""

        # Build label components
        label_parts = []

        if subfamily:
            label_parts.append(subfamily)

        if tribe:
            label_parts.append(tribe)

        label_parts.append(species_name)

        full_label = "_".join(label_parts)

        acc_to_species[acc] = full_label

# ---------- ACCESSION REGEX ----------
acc_pattern = re.compile(r'([A-Z]{1,3}_\d+\.\d+|[A-Z]{2,}\d+\.\d+)')

# ---------- FUNCTION ----------
def relabel_tree(tree_path, output_path):
    tree = Phylo.read(tree_path, "newick")

    for clade in tree.get_terminals():
        label = clade.name
        if not label:
            continue

        acc_match = acc_pattern.search(label)
        if acc_match:
            acc = acc_match.group(1)
            if acc in acc_to_species:
                clade.name = acc_to_species[acc]

    Phylo.write(tree, output_path, "newick")

# ---------- WALK DIRECTORIES ----------
for root, dirs, files in os.walk("."):
    for d in dirs:
        if d.endswith("best_hsp.fasta"):
            dir_path = Path(root) / d

            for file in dir_path.iterdir():
                if file.suffix == ".treefile":

                    input_tree = file
                    output_tree = output_dir / f"labelled_{file.name}"

                    print(f"Processing: {input_tree}")
                    relabel_tree(input_tree, output_tree)

print("\n✅ All done! Labelled trees are in:", output_dir)
