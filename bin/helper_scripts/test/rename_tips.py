#!/usr/bin/env python3


import re
import sys
from Bio import Phylo
sys.setrecursionlimit(100000)

treefile, mapfile = sys.argv[1:3]

mapping = dict(
    line.rstrip("\n").split("\t", 1)
    for line in open(mapfile)
)

tree = Phylo.read(treefile, "newick")

for clade in tree.find_clades():
    if clade.name in mapping:
        clade.name = mapping[clade.name]

Phylo.write(tree, sys.stdout, "newick")
