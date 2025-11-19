# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 11:49:26 2025

@author: lea

python RemoveBirdScaffolds.py birds_scaffolds_file fna_file new_fna_file

"""

import sys
from pathlib import Path


bird_scaffolds_file = Path(sys.argv[1])
fna_file = Path(sys.argv[2])
new_fna_file = Path(sys.argv[3])

bird_scaffolds_list = []

# Read scaffold file and extract contig names
with open(bird_scaffolds_file, "r") as scaffold_file:
    for line in scaffold_file:
        if line.startswith(">"):
            words = line.split()
            contig = words[0][1:]  # Extract contig name after ">"
            bird_scaffolds_list.append(contig)

# Filter the fna file and create a new copy
with open(fna_file, "r") as fna, open(new_fna_file, "w") as newfna:
    keep = True
    for line in fna:
        if line.startswith(">"):
            contig_name = line.split("scaffold=")[-1].split()[0]  # Extract contig name
            keep = contig_name not in bird_scaffolds_list  # Check if it should be removed
        if keep:
            newfna.write(line)
