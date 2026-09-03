#!/usr/bin/env python3
import sys
from ete3 import Tree
from pathlib import Path

if len(sys.argv) != 2:
    print("Usage: python script.py <tree_path>")
    sys.exit(1)

tree_path = sys.argv[1]
path = Path(tree_path)

output_path = path.parent / f"{path.stem}_midcollapse{path.suffix}"

t = Tree(tree_path)
t.set_outgroup(t.get_midpoint_outgroup())

for node in t.traverse():
    if not node.is_leaf() and node.support < 0.5:
        node.delete()

t.write(outfile=str(output_path))
