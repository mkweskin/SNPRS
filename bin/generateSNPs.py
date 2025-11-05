#!/usr/bin/env python3
import pandas as pd
import os
import sys
from ete3 import Tree
from natsort import natsorted
from collections import Counter
import argparse
import polars as pl
import pyarrow.parquet as pq
import shutil
import json

def parse_args():
    parser = argparse.ArgumentParser(description="Create SNPRS SNPs file from tree and grouped data")
    
    # Data args
    parser.add_argument("--out", dest="output_directory", type=str, required=True,help="Path to save output")
    parser.add_argument("--snp_id", dest="snp_id", type=str, required=True,help="Prefix for output files")
    parser.add_argument("--tree", dest="tree_file", type=str, required=True,help="Path to tree file used to generate --splits")
    parser.add_argument("--groups", dest="group_file", type=str, required = True, help="Path to group file created from --tree")
    parser.add_argument("--bases", dest="base_file", type=str, required=True,help="Path to _Bases.parquet")
    parser.add_argument("--scaffold", dest="scaffold_file", type=str, required=True,help="Path to _Scaffold.parquet")

    return parser.parse_args()

def flatten_tips(df):
    return set(
        tip.strip()
        for clade in df["Clade_Tips"].dropna()
        for tip in clade.split(";")
        if tip.strip()
    )

def get_ingroup_parquet(base_file, group_ids, group_name, temp_directory):

    ingroup_parquet = os.path.join(temp_directory, f"{group_name}_Ingroup.parquet")
    if os.path.exists(ingroup_parquet):
        return ingroup_parquet,pq.ParquetFile(ingroup_parquet).metadata.num_rows

    group_ids = list(group_ids)
    group_name = str(group_name)

    lazy_base = pl.scan_parquet(base_file).select(group_ids)
    row_count = pq.ParquetFile(base_file).metadata.num_rows
    lazy_row_numbers = pl.LazyFrame({"row_nr": list(range(row_count))})

    filtered_base_rows = (
        pl.concat([lazy_row_numbers, lazy_base], how="horizontal")
        .filter(pl.all_horizontal(pl.col(group_ids).is_in(["A", "C", "T", "G", "-", "?"])))
        .with_columns([
            pl.concat_list([pl.when(pl.col(c) != "?").then(pl.col(c)).otherwise(None) for c in group_ids])
            .alias("nonzero_values")
        ])
        .with_columns([
            pl.col("nonzero_values").list.drop_nulls().list.unique().alias("unique_values"),
            pl.col("nonzero_values").list.drop_nulls().list.len().alias("In_Count")
        ])
        .filter(pl.col("unique_values").list.len() == 1) # <---- Adjustment point to allow for heterozygous In_Base
        .with_columns([
            pl.col("unique_values").list.first().alias("In_Base")
        ])
        .select(["row_nr", "In_Base", "In_Count"])
    )

    output_rows = filtered_base_rows.select(pl.len()).collect().item()
    if output_rows == 0:
        return None,0

    filtered_base_rows.collect(streaming=True).write_parquet(ingroup_parquet)
    return ingroup_parquet,output_rows

def get_degenerate_outgroup_parquet(ingroup_fixed_parquet,ingroup_id,base_file,group_ids,group_name,temp_directory,higher_level_snp_parquet=None):
    
    degenerate_map = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "-": ["-"],
    "R": ["A","G"],
    "Y": ["C","T"],
    "S": ["G","C"],
    "W": ["A","T"],
    "K": ["G","T"],
    "M": ["A","C"],
    "B": ["C","G","T"],
    "D": ["A","G","T"],
    "H": ["A","C","T"],
    "V": ["A","C","G"],
    }

    degenerate_map.update({k.lower(): v + ["-"] for k, v in degenerate_map.items() if k != "-"})

    fixed_bases = {'A','C','T','G','-','?'}

    outgroup_parquet = os.path.join(temp_directory, f"{ingroup_id}_{group_name}_Outgroup.parquet")

    if os.path.exists(outgroup_parquet):
        return outgroup_parquet,pq.ParquetFile(outgroup_parquet).metadata.num_rows

    group_ids = list(group_ids)
    group_name = str(group_name)

    lazy_ingroup_rows = pl.scan_parquet(ingroup_fixed_parquet).select(['row_nr'])

    lazy_base = pl.scan_parquet(base_file).select(group_ids)
    row_count = pq.ParquetFile(base_file).metadata.num_rows
    lazy_row_numbers = pl.LazyFrame({"row_nr": list(range(row_count))})

    ingroup_base_rows = (
        pl.concat([lazy_row_numbers, lazy_base], how="horizontal")
        .join(lazy_ingroup_rows, on="row_nr", how="inner")
        .filter(~pl.any_horizontal(pl.col(group_ids) == "N"))
    )

    if higher_level_snp_parquet:
        high_frames = [pl.scan_parquet(parquet).select("row_nr") for parquet in higher_level_snp_parquet]
        
        if high_frames:
            lazy_high = pl.concat(high_frames).unique(subset=["row_nr"]).sort("row_nr")
            ingroup_base_rows = ingroup_base_rows.join(lazy_high, on="row_nr", how="anti")

    all_fixed_mask = pl.all_horizontal(pl.col(group_ids).is_in(fixed_bases))
    
    simple_rows = ingroup_base_rows.filter(all_fixed_mask)
    degenerate_rows = ingroup_base_rows.filter(~all_fixed_mask)
    
    simple_outgroup_rows = pl.LazyFrame()
    degenerate_outgroup_rows = pl.LazyFrame()

    simple_outgroup_rows = (
        simple_rows
        .with_columns([
            pl.concat_list([
                pl.when(pl.col(c) != "?").then(pl.col(c)).otherwise(None)
                for c in group_ids
            ]).alias("nonzero_values")
        ])
        .with_columns([
            pl.col("nonzero_values").list.drop_nulls().list.len().alias("Out_Count"),
            pl.col("nonzero_values").list.drop_nulls().list.unique().list.join("").alias("Out_Base"),
        ])
        .filter(pl.col("Out_Count") > 0)
        .select(["row_nr", "Out_Base", "Out_Count"])
    )

    degenerate_outgroup_rows = (
        degenerate_rows
        .with_columns([
            pl.concat_list([
                pl.when(pl.col(c) != "?").then(pl.col(c)).otherwise(None)
                for c in group_ids
            ]).alias("nonzero_values")
        ])
        .with_columns([
            pl.col("nonzero_values").list.drop_nulls().list.len().alias("Out_Count"),
            pl.col("nonzero_values")
            .list.drop_nulls()
            .map_elements(lambda vals: sorted(set(sum((degenerate_map.get(v, []) for v in vals), []))))
            .list.join("")
            .alias("Out_Base"),
        ])
        .filter(pl.col("Out_Count") > 0)
        .select(["row_nr", "Out_Base", "Out_Count"])
    )

    simple_row_count = simple_outgroup_rows.select(pl.len()).collect().item()
    degen_row_count = degenerate_outgroup_rows.select(pl.len()).collect().item()
    total_rows = simple_row_count + degen_row_count

    if total_rows == 0:
        return None, 0

    frames = []
    if simple_row_count > 0:
        frames.append(simple_outgroup_rows)
    if degen_row_count > 0:
        frames.append(degenerate_outgroup_rows)

    (
        pl.concat(frames, how="vertical")
        .sort("row_nr")
        .collect(streaming=True)
        .write_parquet(outgroup_parquet)
    )

    return outgroup_parquet, total_rows

def get_ingroup_snps(ingroup_parquet,outgroup_parquet,ingroup_id,outgroup_id,temp_directory):
        
    lazy_ingroup = pl.scan_parquet(ingroup_parquet)
    lazy_outgroup = pl.scan_parquet(outgroup_parquet)

    ingroup_row_count = pq.ParquetFile(ingroup_parquet).metadata.num_rows
    outgroup_row_count = pq.ParquetFile(outgroup_parquet).metadata.num_rows

    raw_snp_parquet = os.path.join(temp_directory,f"{ingroup_id}_{outgroup_id}_SNPs.parquet")
    
    if os.path.exists(raw_snp_parquet):
        os.remove(raw_snp_parquet)

    joined = (
        lazy_ingroup
        .join(lazy_outgroup, on="row_nr", how="inner")
        .filter(~pl.struct(["In_Base", "Out_Base"])
        .map_elements(lambda x: x["In_Base"] in x["Out_Base"],return_dtype=pl.Boolean))
    )

    joined_row_count = joined.select(pl.len()).collect().item()

    if joined_row_count == 0:
        return None,0

    (
        joined
        .collect(streaming=True)
        .write_parquet(raw_snp_parquet)
    )

    return raw_snp_parquet,joined_row_count

def compileSNPs(non_zero_comparisons, snp_id,output_directory, temp_directory,base_parquet,scaffold_parquet,tree_file,group_file):

    output_file = os.path.join(output_directory,f"{snp_id}_SNPs.parquet")
    output_json = os.path.join(output_directory,f"{snp_id}.json")
    out_csv = os.path.join(output_directory,f"{snp_id}_Comparisons.csv")

    row_count = pq.ParquetFile(scaffold_parquet).metadata.num_rows
    lazy_row_numbers = pl.LazyFrame({"row_nr": list(range(row_count))})
    lazy_scaffold = pl.scan_parquet(scaffold_parquet)
    
    rowed_scaffold = pl.concat([lazy_row_numbers,lazy_scaffold],how="horizontal")

    comparisons = [ (row['Ingroup'],row['Comparison']) for _, row in non_zero_comparisons.iterrows()]
    
    lazy_list = []
    
    for ingroup, outgroup in comparisons:
        
        if outgroup == "No_Outgroup":
            parquet_file = os.path.join(temp_directory, "Ingroup_Ingroup.parquet")
        
        else:
            parquet_file = os.path.join(temp_directory, f"{ingroup}_{outgroup}_SNPs.parquet")

        if not os.path.exists(parquet_file):
            raise ValueError(f"{parquet_file} does not exist...")

        if outgroup == "No_Outgroup":
            lazy_pq = (
                pl.scan_parquet(parquet_file)
                .with_columns([
                    pl.lit(ingroup).alias("Ingroup"),
                    pl.lit(outgroup).alias("Outgroup"),
                    pl.lit("").alias("Out_Base"),
                    pl.col("In_Count").cast(pl.Int32).alias("Out_Count")])
                .select(['row_nr','Ingroup','Outgroup','In_Base','Out_Base','In_Count','Out_Count'])
            ) 
        else:

            lazy_pq = ( 
                pl.scan_parquet(parquet_file)
                .with_columns([
                    pl.lit(ingroup).alias("Ingroup"),
                    pl.lit(outgroup).alias("Outgroup"),
                    pl.col("Out_Count").cast(pl.Int32).alias("Out_Count")])
                .select(['row_nr','Ingroup','Outgroup','In_Base','Out_Base','In_Count','Out_Count'])
            )

        lazy_list.append(lazy_pq)

    if not lazy_list:
        raise ValueError(f"No SNPs could be processed..")

    combined_lazy = pl.concat(lazy_list).sort('row_nr').join(rowed_scaffold,on="row_nr",how="left").sort(['contig_index','contig_position'])

    combined_lazy_count = combined_lazy.select(pl.len()).collect().item()
    combined_unique_count = combined_lazy.select(["contig_index", "contig_position"]).unique().select(pl.len()).collect().item()

    if combined_lazy_count == 0:
        raise ValueError(f"No SNPs could be processed..")
    (
        combined_lazy
        .select(['contig_index','contig_position','Ingroup','Outgroup','In_Base','Out_Base','In_Count','Out_Count'])
        .collect(streaming=True)
        .write_parquet(output_file,compression="snappy")
    )

    json_info = {
        "SNP_ID":snp_id,
        "Base_Parquet":base_parquet,
        "Scaffold_Parquet":scaffold_parquet,
        "Tree_File":tree_file,
        "Groups_File":group_file,
        "SNP_Parquet":output_file,
        "Comparison_Table":out_csv,
        "Total_Rows":str(combined_lazy_count),
        "Unique_Positions":str(combined_unique_count)
        }

    with open(output_json, "w", encoding="utf-8") as f:
        json.dump(json_info, f, indent=4)

    shutil.rmtree(temp_directory)

if __name__ == "__main__":

    args = parse_args()

    # region 00: Parse args

    snp_id = args.snp_id
    output_directory = os.path.abspath(args.output_directory)
    temp_directory = os.path.join(output_directory,"Temp_SNP")

    if os.path.exists(temp_directory):
        shutil.rmtree(temp_directory)
    os.mkdir(temp_directory)

    # Process tree file
    tree_file = os.path.abspath(args.tree_file)
    tree = Tree(tree_file, format=1)
    tree_tips = set(tree.get_leaf_names())

    # Process group file
    group_file = os.path.abspath(args.group_file)
    group_df = pd.read_csv(group_file, sep=",")
    required_cols = {"Clade_ID", "Clade_Tips", "Clade_Type"}
    missing = required_cols - set(group_df.columns)
    if missing:
        raise ValueError(f"Missing columns in input: {', '.join(missing)}")
    group_tips = flatten_tips(group_df)
    assert group_tips == tree_tips, "Group taxa do not match tree taxa"

    duplicates = group_df["Clade_ID"][group_df["Clade_ID"].duplicated()]
    assert duplicates.empty, f"Duplicate Clade_IDs found: {duplicates.tolist()}"

    # Process parquet data
    base_parquet = os.path.abspath(args.base_file)
    scaffold_parquet = os.path.abspath(args.scaffold_file)
    parquet_ids = pl.scan_parquet(base_parquet).collect_schema().names()
    assert set(parquet_ids) == tree_tips, "Parquet taxa do not match tree taxa"

    # endregion

    # region 01: Check outgroup

    outgroup_mask = group_df["Clade_Type"] == "Outgroup"
    ingroup_mask = group_df["Clade_Type"] == "Ingroup"

    outgroup_count = outgroup_mask.sum()
    ingroup_count = ingroup_mask.sum()

    assert outgroup_count <= 1, f"'Outgroup' occurs {outgroup_count} times (expected 0 or 1)."
    assert ingroup_count <= 1, f"'Ingroup' occurs {ingroup_count} times (expected 0 or 1)."

    ingroup_ids = {}
    outgroup_ids = {}
    ingroup_id = "Ingroup"

    if outgroup_count == 1:
        
        outgroup_row = group_df[outgroup_mask].iloc[0]
        outgroup_tip_string = outgroup_row["Clade_Tips"]
        outgroup_ids = {t for t in outgroup_tip_string.split(";")}

        ingroup_row = group_df[ingroup_mask].iloc[0]
        ingroup_tip_string = ingroup_row["Clade_Tips"]
        ingroup_ids = {t for t in ingroup_tip_string.split(";")}
        ingroup_id = ingroup_row["Clade_ID"].strip()

    else:
        ingroup_ids = tree_tips

    if len(ingroup_ids) < 1:
        raise ValueError("No ingroup tips detected")

    # endregion

    # region 02: Process outgroup
    all_comparisons = {}

    ingroup_fixed_parquet, ingroup_fixed_sites = get_ingroup_parquet(base_parquet, ingroup_ids, ingroup_id, temp_directory)

    if ingroup_fixed_sites > 0:
        
        if outgroup_ids:

            overlap = set(ingroup_ids) & set(outgroup_ids)
            if overlap:
                raise ValueError(f"Overlap detected between ingroup and outgroup IDs: {overlap}.")

            assert (set(ingroup_ids) | set(outgroup_ids)) == tree_tips, "Ingroup/Outgroup do not add up to all tips."

            outgroup_degen_parquet, outgroup_count = get_degenerate_outgroup_parquet(ingroup_fixed_parquet, ingroup_id, base_parquet, outgroup_ids, "Outgroup", temp_directory)

            if outgroup_count > 0:
                ingroup_snp_parquet, ingroup_snp_count = get_ingroup_snps(ingroup_fixed_parquet, outgroup_degen_parquet, ingroup_id, "Outgroup", temp_directory)
                all_comparisons[(ingroup_id, "Outgroup")] = (ingroup_fixed_sites,outgroup_count,ingroup_snp_count)
            else:
                all_comparisons[(ingroup_id, "Outgroup")] = (ingroup_fixed_sites,0,0)
        else:
            all_comparisons[(ingroup_id, "No_Outgroup")] = (ingroup_fixed_sites, 0, ingroup_fixed_sites)

    
    else:
        if outgroup_ids:
            all_comparisons[(ingroup_id, "NO_FIXED_SITES")] = (0, 0, 0)
        else:
            all_comparisons[(ingroup_id, "NO_FIXED_SITES")] = (0, 0, 0)

    # endregion

    # region 03: Process ingroups

    ingroup_df = group_df.loc[~group_df["Clade_Type"].isin(["Ingroup", "Outgroup"])].copy()

    if ingroup_df.empty:
        raise ValueError("No ingroup rows remain after processing outgroup.")

    clades = {row["Clade_ID"]: set(row["Clade_Tips"].split(";")) for row in ingroup_df.to_dict(orient="records")}
    
    clade_ids = clades.keys()
    clade_ids_small = sorted(clade_ids, key=lambda c: len(clades[c]))
    clade_ids_large = sorted(clade_ids, key=lambda c: len(clades[c]), reverse=True)

    for clade in clade_ids_small:
        

        ingroup_snps = []

        clade_taxa = set(clades[clade])
        primary_ingroup_fixed_parquet,primary_fixed_count = get_ingroup_parquet(base_parquet, clade_taxa, clade, temp_directory)

        if primary_fixed_count == 0:
            all_comparisons[(clade, "NO_FIXED_SITES")] = (0, 0, 0)
            
        else:
        
            # First get sites that separate the clade from all other ingroup taxa
            outgroup_taxa = set(ingroup_ids) - clade_taxa

            if len(outgroup_taxa) > 0:
                
                outgroup_degen_parquet,outgroup_count = get_degenerate_outgroup_parquet(primary_ingroup_fixed_parquet, clade, base_parquet,outgroup_taxa, ingroup_id, temp_directory)
                
                if outgroup_count > 0:
                    
                    primary_ingroup_snp_parquet,snp_count = get_ingroup_snps(primary_ingroup_fixed_parquet, outgroup_degen_parquet, clade, ingroup_id, temp_directory)

                    all_comparisons[(clade, "Ingroup")] = (primary_fixed_count,outgroup_count,snp_count)
                    
                    if snp_count > 0:
                        ingroup_snps.append(primary_ingroup_snp_parquet)

                else:
                    all_comparisons[(clade, "Ingroup")] = (primary_fixed_count,0,0)

            # Inner loop also needs descending order
            for clade2 in clade_ids_large:
                
                if clade == clade2:
                    continue

                clade2_taxa = set(clades[clade2])

                if clade2_taxa == ingroup_ids or clade2_taxa == outgroup_taxa:
                    continue

                if clade_taxa < clade2_taxa:
            
                    sub_outgroup_taxa = clade2_taxa - clade_taxa
                    outgroup_degen_parquet,outgroup_count = get_degenerate_outgroup_parquet(primary_ingroup_fixed_parquet, clade, base_parquet,sub_outgroup_taxa, clade2, temp_directory,ingroup_snps)
                    
                    if outgroup_count > 0:
                        
                        ingroup_snp_parquet,snp_count = get_ingroup_snps(primary_ingroup_fixed_parquet, outgroup_degen_parquet, clade, clade2, temp_directory)

                        all_comparisons[(clade, clade2)] = (primary_fixed_count,outgroup_count,snp_count)
                        
                        if snp_count > 0:
                            ingroup_snps.append(ingroup_snp_parquet)

                    else:
                        all_comparisons[(clade, clade2)] = (primary_fixed_count,0,0)


    all_comparisons_df = pd.DataFrame(
        [
            {
                "Ingroup": k[0],
                "Comparison": k[1],
                "Fixed_Sites": v[0],
                "Outgroup_Count": v[1],
                "SNP_Count": v[2],
            }
            for k, v in all_comparisons.items()
        ]
    )

    all_comparisons_df["Min_Ingroup"] = ""
    all_comparisons_df["Min_Outgroup"] = ""
    all_comparisons_df["Fixed_Outgroup"] = ""

    out_csv = os.path.join(output_directory,f"{snp_id}_Comparisons.csv")
    all_comparisons_df.to_csv(out_csv, sep=",", index=False)
    
    non_zero_comparisons = all_comparisons_df.loc[all_comparisons_df["SNP_Count"] > 0].copy()

    if not non_zero_comparisons.empty:
        compileSNPs(non_zero_comparisons, snp_id, output_directory, temp_directory,base_parquet,scaffold_parquet,tree_file,group_file)
    
    # endregion
