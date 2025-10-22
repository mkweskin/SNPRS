#!/usr/bin/env python3
import pandas as pd
import os
import sys
from ete3 import Tree
from natsort import natsorted
from collections import Counter
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Create SNPRS grouping file from tree and splits data")
    
    # Data args
    parser.add_argument("--tree", dest="tree_file", type=str, required=True,help="Path to tree file used to generate --splits")
    parser.add_argument("--splits", dest="splits_file", type=str, default = None, help="Path to splits file created from --tree, with or without names [Generate from tree if not provided]")
    parser.add_argument("--mono", dest="is_monotypic",action = "store_true",help="Process singleton data from splits file as groups")

    return parser.parse_args()

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
                    "Clade_ID":"",
                    "Mirror_ID":""
                })

    for leaf in tree.iter_leaves():
        rows.append({
            "Clade_Taxa": leaf.name,
            "Mirror_Taxa": "",
            "Clade_ID":"",
            "Mirror_ID":""
        })

    return pd.DataFrame(rows)

def flatten_tips(df):
    return set(
        tip.strip()
        for clade in df["Clade_Tips"].dropna()
        for tip in clade.split(";")
        if tip.strip()
    )

def is_monophyletic(tree, tips):
    try:
        mrca = tree.get_common_ancestor(tips)
    except Exception:
        return False 
    leaf_names = set(leaf.name for leaf in mrca.get_leaves())
    return leaf_names == set(tips)

def split_ingroup_outgroup(split_df,all_taxa,is_blank,only_outgroup,is_mono,tree):

    outgroup_df = pd.DataFrame(columns=["Clade_ID", "Clade_Tips", "Group_Type"])

    for col in ["Clade_ID", "Mirror_ID"]:
        split_df[col] = split_df[col].fillna("").astype(str).str.strip()

    outgroup_mask = (
        (split_df["Clade_ID"].str.lower() == "outgroup") |
        (split_df["Mirror_ID"].str.lower() == "outgroup")
    )

    outgroup_count = outgroup_mask.sum()
    assert outgroup_count <= 1, f"'Outgroup' occurs {outgroup_count} times (expected 0 or 1)."

    ingroup_tips = {}
    outgroup_tips = {}
    ingroup_tip_string = ""
    outgroup_tip_string = ""

    if outgroup_count == 1:
        outgroup_rows = []
        outgroup_row = split_df[outgroup_mask].iloc[0]

        if outgroup_row["Clade_ID"].strip().lower() == "outgroup":
            outgroup_tip_string = outgroup_row["Clade_Taxa"]
            ingroup_id = outgroup_row["Mirror_ID"].strip() or "Ingroup"
        else:
            outgroup_tip_string = outgroup_row["Mirror_Taxa"]
            ingroup_id = outgroup_row["Clade_ID"].strip() or "Ingroup"

        outgroup_tips = {t for t in outgroup_tip_string.split(";")}
        ingroup_tips = {t for t in all_taxa - outgroup_tips}
        
        if len(ingroup_tips) < 1:
            raise ValueError("No ingroup tips detected")

        ingroup_tip_string = ";".join(natsorted(ingroup_tips))
        
        outgroup_rows.append({
            "Clade_ID": ingroup_id,
            "Clade_Tips": ingroup_tip_string,
            "Group_Type": "Ingroup"
        })

        outgroup_rows.append({
            "Clade_ID": "Outgroup",
            "Clade_Tips": outgroup_tip_string,
            "Group_Type": "Outgroup"
        })

        outgroup_df = pd.DataFrame(outgroup_rows)
    else:
        ingroup_tips = all_taxa
    
    if is_mono:
        ingroup_df = (
            split_df
            .loc[
                (split_df["Mirror_Taxa"].astype(str).str.strip() == "") &
                (split_df["Clade_Taxa"].isin(ingroup_tips))
            ]
            .rename(columns={"Clade_Taxa": "Clade_Tips"})
            .assign(
                Clade_ID=lambda df: df["Clade_ID"].astype(str).str.strip().replace("", df["Clade_Tips"])
            )
            [["Clade_ID", "Clade_Tips"]]
        )

        ingroup_df['Group_Type'] = "Mono"

    else:

        ingroup_rows = []
        seen_tips = {str(ingroup_tip_string).strip(), str(outgroup_tip_string).strip()}

        if is_blank or only_outgroup:
            
            id_type = "Blank" if is_blank else "Only_Outgroup"
            
            ingroup_df = (
                split_df
                .loc[
                    (split_df["Mirror_Taxa"].astype(str).str.strip() == "") &
                    (split_df["Clade_Taxa"].isin(ingroup_tips))
                ]
                .rename(columns={"Clade_Taxa": "Clade_Tips"})
                .assign(
                    Clade_ID=lambda df: df["Clade_ID"].astype(str).str.strip().replace("", df["Clade_Tips"]),
                    Group_Type=id_type)
                [["Clade_ID", "Clade_Tips","Group_Type"]]
            )

            all_splits = (
                split_df["Clade_Taxa"].astype(str).str.strip().tolist() +
                split_df["Mirror_Taxa"].astype(str).str.strip().tolist()
            )

            nonempty_splits = [x for x in all_splits if x and ";" in x]
            ingroup_splits = [
                x for x in nonempty_splits
                if all(tip.strip() in ingroup_tips for tip in x.split(";") if tip.strip())
            ]

            i = 1
            for id_col, tips_col in [("Clade_ID", "Clade_Taxa"), ("Mirror_ID", "Mirror_Taxa")]:
                mask = split_df[tips_col].astype(str).str.strip().isin(ingroup_splits)
                for _, row in split_df.loc[mask].iterrows():
                    
                    string_a = str(row[tips_col]).strip()
                    a_tips = {t.strip() for t in string_a.split(";") if t.strip()}
                    
                    b_tips = set(ingroup_tips) - set(a_tips)
                    string_b = ";".join(natsorted(b_tips))

                    if string_a and string_a not in seen_tips:                        
                        ingroup_rows.append({
                            "Clade_ID": f"Auto_Group_{i}",
                            "Clade_Tips": string_a,
                            "Group_Type":id_type
                        })
                        
                        seen_tips.add(string_a)
                        i += 1                    
                    
                    if string_b and string_b not in seen_tips:
                        ingroup_rows.append({
                            "Clade_ID": f"Auto_Group_{i}",
                            "Clade_Tips": string_b,
                            "Group_Type":id_type
                        })
                        seen_tips.add(string_b)
                        i += 1

            auto_rows = pd.DataFrame(ingroup_rows)
            ingroup_df = pd.concat([ingroup_df, auto_rows], ignore_index=True)


        # If more than just outgroup is labeled, all terminal groups must be labeled
        else:

            for id_col, tips_col in [("Clade_ID", "Clade_Taxa"), ("Mirror_ID", "Mirror_Taxa")]:
                mask = (split_df[id_col].str.lower() != "outgroup") & (split_df[id_col] != "")
                for _, row in split_df[mask].iterrows():
                    clade_tips = row[tips_col]
                    if clade_tips and clade_tips not in seen_tips:
                        ingroup_rows.append({
                            "Clade_ID": row[id_col],
                            "Clade_Tips": clade_tips
                        })

        ingroup_df = pd.DataFrame(ingroup_rows)

    return ingroup_df, outgroup_df, ingroup_tips, outgroup_tips

def assign_terminal(ingroup_df, ingroup_ids):
    remaining_ids = set(ingroup_ids)
    terminal_rows = []

    tip_to_rows = {}
    for _, row in ingroup_df.iterrows():
        tips = set(t.strip() for t in row["Clade_Tips"].split(";") if t.strip())
        for tip in tips:
            tip_to_rows.setdefault(tip, []).append((row["Clade_ID"], tips))

    while remaining_ids:
        tip = next(iter(remaining_ids))

        if tip not in tip_to_rows:
            raise ValueError(f"Cannot assign tip '{tip}' to any terminal group. Check groupings/singletons...")
    
        candidate_rows = tip_to_rows[tip]
        lengths = [len(tips) for _, tips in candidate_rows]
        min_len = min(lengths)

        if lengths.count(min_len) > 1:
            tied_ids = [cid for (cid, tips), l in zip(candidate_rows, lengths) if l == min_len]
            raise ValueError(
                f"Tie detected for tip '{tip}': multiple Clade_IDs have minimal size {min_len} "
                f"({tied_ids}). Cannot determine unique terminal group.")

        terminal_clade_id, terminal_tips = min(candidate_rows, key=lambda x: len(x[1]))

        terminal_rows.append({
            "Clade_ID": terminal_clade_id,
            "Clade_Tips": ";".join(natsorted(terminal_tips)),
            "Group_Type": "Terminal"
        })

        remaining_ids -= terminal_tips

    if not terminal_rows:
        raise ValueError(f"Cannot process terminal groups. Check groupings/singletons...")

    terminal_group_df = pd.DataFrame(terminal_rows)

    terminal_tips_flat = [
    tip.strip()
    for clade in terminal_group_df["Clade_Tips"]
    for tip in clade.split(";")
    if tip.strip()
    ]

    terminal_tips_set = set(terminal_tips_flat)

    missing_tips = ingroup_ids - terminal_tips_set
    if missing_tips:
        raise ValueError(f"Some ingroup tips are missing from terminal groups: {missing_tips}")


    tip_counts = Counter(terminal_tips_flat)
    duplicates = [tip for tip, count in tip_counts.items() if count > 1]
    if duplicates:
        raise ValueError(f"Some tips appear more than once in terminal groups: {duplicates}")

    used_ids = set(terminal_group_df["Clade_ID"])
    named_internal_rows = []

    for _, row in ingroup_df.iterrows():
        clade_id = str(row.get("Clade_ID", "")).strip()
        clade_tips = row["Clade_Tips"]

        if clade_id and clade_id not in used_ids:
            named_internal_rows.append({
                "Clade_ID": clade_id,
                "Clade_Tips": clade_tips,
                "Group_Type": "Internal"
            })

    named_internal_df = pd.DataFrame(named_internal_rows)
    return pd.concat([terminal_group_df, named_internal_df])

def add_internals(terminal_named_df, split_df,tree):

    existing_clade_tips = set(
        ";".join(natsorted(t.strip() for t in clade.split(";") if t.strip()))
        for clade in terminal_named_df["Clade_Tips"]
    )

    new_internal_rows = []
    i = 1

    for id_col, tip_col in [("Clade_ID", "Clade_Taxa"), ("Mirror_ID", "Mirror_Taxa")]:
        mask = split_df[id_col].str.lower() != "outgroup"
        for _, row in split_df[mask].iterrows():
            if pd.isna(row[tip_col]) or str(row[tip_col]).strip() == "":
                continue

            split_tips = [t.strip() for t in row[tip_col].split(";") if t.strip()]
            normalized_tips = ";".join(natsorted(split_tips))

            if (normalized_tips not in existing_clade_tips) & (is_monophyletic(tree,split_tips)):
                new_internal_rows.append({
                    "Clade_ID": f"SNPRS_Internal_{i}",
                    "Clade_Tips": normalized_tips,
                    "Group_Type": "Internal"
                })
                existing_clade_tips.add(normalized_tips)
                i += 1

    if not new_internal_rows:
        return terminal_named_df
    else:
        unnamed_internal_df = pd.DataFrame(new_internal_rows)
        terminal_named_df = pd.concat([terminal_named_df, unnamed_internal_df], ignore_index=True)

    return terminal_named_df


if __name__ == "__main__":

    args = parse_args()

    # Process tree file
    tree_file = os.path.abspath(args.tree_file)
    out_csv = os.path.splitext(tree_file)[0] + "_Monophyletic_Groups.csv"
    tree = Tree(tree_file, format=1)
    tree_tips = set(tree.get_leaf_names())

    # Process splits file
    splits_file = os.path.abspath(args.splits_file)

    if not splits_file:
        split_df = extract_splits(tree)
    else:
        split_df = pd.read_csv(splits_file, sep=",")
        required_cols = {"Clade_Taxa", "Mirror_Taxa", "Clade_ID", "Mirror_ID"}
        missing = required_cols - set(split_df.columns)
        if missing:
            raise ValueError(f"Missing columns in input: {', '.join(missing)}")

    singleton_rows = split_df["Mirror_Taxa"].isna() | (split_df["Mirror_Taxa"].astype(str).str.strip() == "")
    csv_taxa = set(split_df[singleton_rows]['Clade_Taxa'])
    assert csv_taxa == tree_tips, "CSV taxa do not match tree taxa"
    
    is_blank = (
        not split_df["Clade_ID"].fillna("").astype(str).str.strip().any() and
        not split_df["Mirror_ID"].fillna("").astype(str).str.strip().any()
    )

    is_mono = args.is_monotypic
    
    # Check if only the outgroup is noted
    all_ids = (
        split_df["Clade_ID"].fillna("").astype(str).str.strip().tolist() +
        split_df["Mirror_ID"].fillna("").astype(str).str.strip().tolist()
    )

    nonempty_ids = [x for x in all_ids if x]

    only_outgroup = (
        len(set(x.lower() for x in nonempty_ids)) == 1 and
        set(x.lower() for x in nonempty_ids) == {"outgroup"}
    )

    ingroup_df,outgroup_df,ingroup_ids,outgroup_ids = split_ingroup_outgroup(split_df,tree_tips,is_blank,only_outgroup,is_mono,tree)
    
    if outgroup_ids:
        overlap = set(ingroup_ids) & set(outgroup_ids)
        if overlap:
            raise ValueError(
                f"Overlap detected between ingroup and outgroup IDs: {overlap}. "
                "Each tip should be assigned to only one group."
            )
    
    assert set(set(ingroup_ids) | set(outgroup_ids)) == tree_tips, "Ingroup/Outgroup to not add up to all tips"

    if is_blank or only_outgroup or is_mono:
        pd.concat([ingroup_df,outgroup_df]).to_csv(out_csv, sep=",", index=False)

    else:
        terminal_named_df = assign_terminal(ingroup_df, ingroup_ids)
        full_ingroup_df = add_internals(terminal_named_df, split_df,tree)
        pd.concat([full_ingroup_df,outgroup_df]).to_csv(out_csv, sep=",", index=False)
