from ete3 import Tree
import sys
import os

tree_file = os.path.abspath(sys.argv[1])

if not os.path.exists(tree_file):
    raise FileNotFoundError(f"Tree file not found: {tree_file}")

try:
    tree = Tree(tree_file, format=1)
except Exception as e:
    raise ValueError(f"Failed to parse tree '{tree_file}': {e}")