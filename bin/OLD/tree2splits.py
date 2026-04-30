#!/usr/bin/env python3
import sys
import os
import pandas as pd
from ete3 import Tree
from natsort import natsorted

def extract_splits(tree):
    all_taxa = set(tree.get_leaf_names())
    rows = []

    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            clade_taxa = set(node.get_leaf_names())

            if clade_taxa != all_taxa:
                mirror_taxa = all_taxa - clade_taxa
                rows.append({
                    "Clade_Taxa": ";".join(natsorted(clade_taxa)),
                    "Mirror_Taxa": ";".join(natsorted(mirror_taxa)),
                })

    for leaf in tree.iter_leaves():
        rows.append({
            "Clade_Taxa": leaf.name,
            "Mirror_Taxa": "",
        })

    return pd.DataFrame(rows)

if __name__ == "__main__":

    treefile = os.path.abspath(sys.argv[1])
    outfile = os.path.abspath(sys.argv[2])

    if not os.path.exists(treefile):
        sys.exit(f"Error: {treefile} not found.")

    tree = Tree(treefile, format=1)

    df = extract_splits(tree)
    df["Clade_ID"] = ""
    df["Mirror_ID"] = ""

    df.to_csv(outfile, sep=",", index=False)
