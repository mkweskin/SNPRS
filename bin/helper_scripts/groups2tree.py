import os
import sys
import pandas as pd
from ete3 import Tree,TreeNode

# Load CSV
if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} <input_csv>")
    sys.exit(1)

input_csv = os.path.abspath(sys.argv[1])
internal_df = pd.read_csv(input_csv)
output_file = os.path.join(os.path.dirname(input_csv), "SNP_Groups.nwk")

# Create leaf nodes for all tips
all_tips = set()
for tips in internal_df['Clade_Tips']:
    all_tips.update(tips.split(';'))

nodes = {tip: TreeNode(name=tip) for tip in all_tips}

# Sort internal nodes by number of tips (smallest first)
internal_df['num_tips'] = internal_df['Clade_Tips'].apply(lambda x: len(x.split(';')))
internal_df = internal_df.sort_values('num_tips')

# Build the tree
for _, row in internal_df.iterrows():
    tips = set(row['Clade_Tips'].split(';'))
    clade_id = row['Clade_ID']

    # Find children: nodes whose leaves are fully contained in this clade
    children = []
    for n_name, n_node in list(nodes.items()):
        node_tips = set(n_node.get_leaf_names())
        if node_tips.issubset(tips):
            children.append(n_node)
            nodes.pop(n_name)  # remove from dict; now attached

    # Create internal node and add children
    internal_node = TreeNode(name=clade_id)
    for child in children:
        internal_node.add_child(child)

    # Add internal node to dict
    nodes[clade_id] = internal_node

# There should be a single node left: the root
roots = list(nodes.values())
if len(roots) > 1:
    root_node = TreeNode(name="Root")
    for r in roots:
        root_node.add_child(r)
else:
    root_node = roots[0]

if 'Group_Type' in internal_df.columns and any(internal_df['Group_Type'].str.lower() == 'outgroup'):
    outgroup_row = internal_df[internal_df['Group_Type'].str.lower() == 'outgroup'].iloc[0]
    outgroup_tips = set(outgroup_row['Clade_Tips'].split(';'))

    # Convert to Tree for rooting support
    t = Tree(root_node.write(format=9))  # format=9 preserves structure & names
    outgroup_leaves = [leaf for leaf in t if leaf.name in outgroup_tips]

    if outgroup_leaves:
        common_ancestor = t.get_common_ancestor(outgroup_leaves)
        t.set_outgroup(common_ancestor)
    else:
        print("Warning: Outgroup tips not found in tree — skipping rooting.")
else:
    t = Tree(root_node.write(format=9))

for _, row in internal_df.iterrows():
    tips = set(row['Clade_Tips'].split(';'))
    clade_id = row['Clade_ID']

    for node in t.traverse("postorder"):
        node_tips = set(node.get_leaf_names())
        if node_tips == tips:
            node.name = clade_id
            break

with open(output_file, "w") as f:
    f.write(t.write(format=1) + "\n")