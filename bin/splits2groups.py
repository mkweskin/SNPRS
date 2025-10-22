#!/usr/bin/env python3
import pandas as pd
import os
import sys
from ete3 import Tree
from natsort import natsorted
from collections import Counter


def is_monophyletic(tree, tips):
    try:
        mrca = tree.get_common_ancestor(tips)
    except Exception:
        return False 
    leaf_names = set(leaf.name for leaf in mrca.get_leaves())
    return leaf_names == set(tips)

def split_ingroup_outgroup(split_df):

    for col in ["Clade_ID", "Mirror_ID"]:
        split_df[col] = split_df[col].fillna("").astype(str).str.strip()

    outgroup_mask = (
        (split_df["Clade_ID"].str.lower() == "outgroup") |
        (split_df["Mirror_ID"].str.lower() == "outgroup")
    )

    outgroup_count = outgroup_mask.sum()
    assert outgroup_count <= 1, f"'Outgroup' occurs {outgroup_count} times (expected 0 or 1)."

    outgroup_df = pd.DataFrame(columns=["Clade_ID", "Clade_Tips", "Group_Type"])
    ingroup_tips = ""
    outgroup_tips = ""

    if outgroup_count == 1:
        outgroup_rows = []
        outgroup_row = split_df[outgroup_mask].iloc[0]

        if outgroup_row["Clade_ID"].strip().lower() == "outgroup":
            outgroup_tips = outgroup_row["Clade_Taxa"]
            ingroup_tips = outgroup_row["Mirror_Taxa"]
            ingroup_id = outgroup_row["Mirror_ID"].strip() or "Ingroup"
        else:
            outgroup_tips = outgroup_row["Mirror_Taxa"]
            ingroup_tips = outgroup_row["Clade_Taxa"]
            ingroup_id = outgroup_row["Clade_ID"].strip() or "Ingroup"

        outgroup_rows.append({
            "Clade_ID": ingroup_id,
            "Clade_Tips": ingroup_tips,
            "Group_Type": "Ingroup"
        })

        outgroup_rows.append({
            "Clade_ID": "Outgroup",
            "Clade_Tips": outgroup_tips,
            "Group_Type": "Outgroup"
        })

        outgroup_df = pd.DataFrame(outgroup_rows)

    seen_tips = {ingroup_tips, outgroup_tips}
    ingroup_rows = []

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

    def flatten_tips(df):
        return set(
            tip.strip()
            for clade in df["Clade_Tips"].dropna()
            for tip in clade.split(";")
            if tip.strip()
        )

    ingroup_ids = flatten_tips(ingroup_df)
    outgroup_ids = flatten_tips(outgroup_df[outgroup_df["Group_Type"] == "Outgroup"]) if not outgroup_df.empty else set()

    return ingroup_df, outgroup_df, ingroup_ids, outgroup_ids

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

    new_sets = []
    for id_col, tip_col in [("Clade_ID", "Clade_Taxa"), ("Mirror_ID", "Mirror_Taxa")]:
        mask = split_df[id_col].str.lower() != "outgroup"
        for _, row in split_df[mask].iterrows():
            if pd.isna(row[tip_col]) or str(row[tip_col]).strip() == "":                
                continue
            split_tips = natsorted(t.strip() for t in row[tip_col].split(";") if t.strip())
            normalized_tips = ";".join(split_tips)
            if (normalized_tips not in existing_clade_tips) & (is_monophyletic(tree,split_tips)):
                new_sets.append({
                    "Clade_ID": row[id_col],
                    "Clade_Tips": normalized_tips
                })

    if not new_sets:
        return terminal_named_df

    new_df = pd.DataFrame(new_sets)

    new_internal_rows = []
    i = 1

    for id_col, tip_col in [("Clade_ID", "Clade_Taxa"), ("Mirror_ID", "Mirror_Taxa")]:
        mask = split_df[id_col].str.lower() != "outgroup"
        for _, row in split_df[mask].iterrows():
            if pd.isna(row[tip_col]) or str(row[tip_col]).strip() == "":
                continue

            split_tips = [t.strip() for t in row[tip_col].split(";") if t.strip()]
            normalized_tips = ";".join(natsorted(split_tips))

            if normalized_tips not in existing_clade_tips and is_monophyletic(tree, split_tips):
                new_internal_rows.append({
                    "Clade_ID": f"SNPRS_Internal_{i}",
                    "Clade_Tips": normalized_tips,
                    "Group_Type": "Internal"
                })
                existing_clade_tips.add(normalized_tips)
                i += 1

    if new_internal_rows:
        unnamed_internal_df = pd.DataFrame(new_internal_rows)
        terminal_named_df = pd.concat([terminal_named_df, unnamed_internal_df], ignore_index=True)

    return terminal_named_df


if __name__ == "__main__":

    splits_file = os.path.abspath(sys.argv[1])
    tree_file = os.path.abspath(sys.argv[2])
    tree = Tree(tree_file, format=1)
    all_taxa = set(tree.get_leaf_names())

    split_df = pd.read_csv(splits_file, sep=",")
    required_cols = {"Clade_Taxa", "Mirror_Taxa", "Clade_ID", "Mirror_ID"}
    missing = required_cols - set(split_df.columns)
    if missing:
        raise ValueError(f"Missing columns in input: {', '.join(missing)}")

    singleton_rows = split_df["Mirror_Taxa"].isna() | (split_df["Mirror_Taxa"].astype(str).str.strip() == "")
    csv_taxa = set(split_df[singleton_rows]['Clade_Taxa'])
    assert csv_taxa == all_taxa, "CSV taxa do not match tree taxa"

    ingroup_df,outgroup_df,ingroup_ids,outgroup_ids = split_ingroup_outgroup(split_df)
    
    if outgroup_ids:
        overlap = ingroup_ids & outgroup_ids
        if overlap:
            raise ValueError(
                f"Overlap detected between ingroup and outgroup IDs: {overlap}. "
                "Each tip should be assigned to only one group."
            )
    
    assert set(ingroup_ids | outgroup_ids) == all_taxa, "Ingroup/Outgroup to not add up to all tips"

    terminal_named_df = assign_terminal(ingroup_df, ingroup_ids)
    full_ingroup_df = add_internals(terminal_named_df, split_df,tree)
    output_df = pd.concat([full_ingroup_df,outgroup_df])
    
    out_csv = os.path.splitext(tree_file)[0] + "_Monophyletic_Groups.csv"
    output_df.to_csv(out_csv, sep=",", index=False)
