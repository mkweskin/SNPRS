#!/usr/bin/env python3

import polars as pl
import os
import sys
import argparse
import json
from Bio import SeqIO,AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from natsort import natsorted
import pandas as pd
import math



def getFixedSites(code_parquet, sample_list, group_name):
    fixed_set = {1,2,3,4,16,33,34,35,36,48}

    lazy = pl.scan_parquet(code_parquet).select(
        ['contig_index', 'contig_position'] + sample_list
    )

    all_missing = (
        lazy
        .filter(pl.all_horizontal([pl.col(s) == 0 for s in sample_list]))
        .with_columns([
            pl.lit(0).cast(pl.Int32).alias(group_name),
            pl.lit(len(sample_list)).cast(pl.Int32).alias('Fixed_Count')
        ])
        .select([
            "contig_index",
            "contig_position",
            group_name,
            "Fixed_Count"
        ])
    )

    fixed = (
        lazy
        .filter(~pl.all_horizontal([pl.col(s) == 0 for s in sample_list]))
        .filter(
            pl.max_horizontal([
                pl.when(pl.col(s) != 0).then(pl.col(s)).otherwise(None) 
                for s in sample_list
            ]) ==
            pl.min_horizontal([
                pl.when(pl.col(s) != 0).then(pl.col(s)).otherwise(None) 
                for s in sample_list
            ])
        )
        .with_columns([
            pl.max_horizontal([pl.col(s) for s in sample_list]).cast(pl.Int32).alias(group_name)
        ])
        .filter(pl.col(group_name).is_in(fixed_set))
        .with_columns([
            pl.sum_horizontal([(pl.col(s) != 0).cast(pl.Int32) for s in sample_list]).alias("Fixed_Count")
        ])
        .select(['contig_index','contig_position', group_name, 'Fixed_Count'])
    )


    result = (
        pl.concat([all_missing, fixed])
        .select([
            "contig_index",
            "contig_position",
            group_name,
            "Fixed_Count"
        ])
        .sort(['contig_index','contig_position'])
        .collect(engine="streaming")
    )

    return result

#####

def snpSubtractor(fixed_df, focal_id, code_file, subtract_samples):

    fixed_codes = {1, 2, 3, 4, 16}

    degenerate_map = {
        1:  {5,21,6,22,7,23,11,27,12,28,13,29,15,17,31},
        2:  {5,21,8,24,9,25,11,27,12,28,14,30,15,18,31},
        3:  {6,22,8,24,10,26,11,27,13,29,14,30,15,19,31},
        4:  {7,23,9,25,10,26,12,28,13,29,14,30,15,20,31},
        16: {17,18,29,20,21,22,23,24,25,26,27,28,29,30,31},
    }

    def normalize_code(col):
        return pl.when(col >= 33).then(col - 32).otherwise(col)

    def degenerate_match_expr(focal_col, sample_col):
        return (
            pl.when(focal_col == 1).then(sample_col.is_in(degenerate_map[1]))
            .when(focal_col == 2).then(sample_col.is_in(degenerate_map[2]))
            .when(focal_col == 3).then(sample_col.is_in(degenerate_map[3]))
            .when(focal_col == 4).then(sample_col.is_in(degenerate_map[4]))
            .when(focal_col == 16).then(sample_col.is_in(degenerate_map[16]))
            .otherwise(False)
        )

    focal_df = (
        fixed_df
        .filter(pl.col(focal_id) > 0)
        .with_columns(normalize_code(pl.col(focal_id)).alias(focal_id))
    )
    
    site_count = focal_df.height
    lazy_coords = pl.DataFrame(focal_df.select(["contig_index", "contig_position"])).lazy()

    lazy_code = pl.scan_parquet(code_file)

    for sample in subtract_samples:

        if site_count == 0:
            break

        sample_df = (
            lazy_code
            .select(["contig_index", "contig_position", sample])
            .join(lazy_coords,on=["contig_index","contig_position"],how="inner")
            .with_columns(normalize_code(pl.col(sample)).alias(sample))
            .collect(engine="streaming")
        )

        compare = focal_df.join(
            sample_df,
            on=["contig_index", "contig_position"],
            how="left"
        )

        s = pl.col(sample)
        f = pl.col(focal_id)

        match_mask = (
            (s.is_in(fixed_codes) & (s == f)) |
            ((~s.is_in(fixed_codes)) & (s > 0) & degenerate_match_expr(f, s)) |
            (s < 0)
        )

        matched = compare.filter(match_mask).select(
            ["contig_index", "contig_position"]
        )

        focal_df = focal_df.join(matched, on=["contig_index","contig_position"], how="anti")

        site_count = focal_df.height
        lazy_coords = pl.DataFrame(focal_df.select(["contig_index", "contig_position"])).lazy()


    return focal_df
    """
    for sample in subtract_samples:

        if site_count == 0:
            break

        sample_df = (
            lazy_code
            .select(["contig_index", "contig_position", sample])
            .with_columns(normalize_code(pl.col(sample)).alias(sample))
        ).collect(engine="streaming")

        compare_df = focal_df.join(
            sample_df, on=["contig_index", "contig_position"], how="left"
        )

        sample_col = pl.col(sample)
        focal_col  = pl.col(focal_id)

        fixed_df2 = compare_df.filter(sample_col.is_in(fixed_codes)).with_columns(
            (focal_col == sample_col).alias("Match")
        )

        het_df2 = compare_df.filter((~sample_col.is_in(fixed_codes)) & (sample_col > 0)).with_columns(
            degenerate_match_expr(focal_col, sample_col).alias("Match")
        )

        ploidy_df = compare_df.filter(sample_col < 0)

        match_df = pl.concat([
            fixed_df2.filter(pl.col("Match")).select(["contig_index", "contig_position"]),
            het_df2.filter(pl.col("Match")).select(["contig_index", "contig_position"]),
            ploidy_df.select(["contig_index", "contig_position"])
        ])

        focal_df = focal_df.join(match_df, on=["contig_index","contig_position"], how="anti")

    return focal_df
"""
#####

parser = argparse.ArgumentParser(description='Generate alignment from SNPRS data')
parser.add_argument('-j','--json_file',dest="json_file", type=str,required=True, help='Path to SNPRS joined JSON file')
parser.add_argument('-o','--out',dest="out_dir", type=str,required=True, help='Path to output directory')
parser.add_argument('-n','--name',dest="snp_name", required=True,type=str, help='Prefix for output files')
parser.add_argument('-g','--groups',dest="group_file", required = True, type=str, help='Path to TSV file with group information (Sample_ID, SNP_Group)')

parser.add_argument('-p','--prop',dest="prop", default=None,type=float, help='Proportion of samples within a group required to call a SNP [Default: Use all sites]')
parser.add_argument('-t','--top',dest="top", default=None,type=int, help='Choose the top <top> SNPs based on total missing data')
parser.add_argument('--min_table',dest="min_table", default=None,type=str, help='Path to TSV file with group count information (SNP_Group,Min_Samples)')
parser.add_argument('--missing',dest="missing", default=None,type=int, help='Max number of samples with missing data overall to call a SNP [Default: Use all sites]')

args = parser.parse_args()

json_file = os.path.abspath(args.json_file)

with open(json_file, "r") as f:
    data = json.load(f)

join_id = data["Join_ID"]
joined_directory = data["Joined_Directory"]
sample_ids = natsorted(data["Sample_IDs"].split(","))
scaffold_file = data["Scaffold_File"]
code_file = data["Code_File"]
site_file = data["Site_File"]
sample_summary_file = data["Sample_Summary_File"]
site_count_file = data["Site_Count_File"]

lazy_codes = pl.scan_parquet(code_file)
lazy_sites = pl.scan_parquet(site_file)
missing_info = lazy_sites.select(['contig_index','contig_position','Missing']).collect(engine="streaming")

# Set sample groups
group_file = os.path.abspath(args.group_file)
group_df = pd.read_csv(group_file,sep="\t")

group_data =  (
    group_df.groupby("SNP_Group")["Sample_ID"]
    .apply(list)
    .to_dict()
)

all_ids = list({id for ids in group_data.values() for id in ids})

# Get fixed sites for each group
output_directory = os.path.abspath(args.out_dir)
fixed_sites_dir = os.path.join(output_directory,"Fixed_Sites")

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

if not os.path.exists(fixed_sites_dir):
    os.mkdir(fixed_sites_dir)
else:
    sys.exit(f"{fixed_sites_dir} already exists...")
    
fixed_sites = [
    getFixedSites(code_file, ids, group)
    for group, ids in group_data.items()
]

for i,(group, ids) in enumerate(group_data.items()):
    out_path = os.path.join(fixed_sites_dir, f"{group}.parquet")
    fixed_sites[i].write_parquet(out_path,compression="snappy")
    

snp_results = []

for i,(group, ids) in enumerate(group_data.items()):
    print(f"Starting {group}...\n")
    non_focal = list(set(all_ids) - set(ids))
    sp_snp_df = snpSubtractor(fixed_sites[i],group,code_file,non_focal)

    snp_results.append(sp_snp_df.with_columns(pl.lit(group).alias("SNP_Group"),
    pl.col(group).alias("SNP_Base")).select(['contig_index','contig_position','SNP_Group','SNP_Base','Fixed_Count']))

full_snp_df = pl.concat(snp_results).join(missing_info,on=["contig_index","contig_position"],how="left")

full_snp_count_df = full_snp_df.group_by("SNP_Group").agg(pl.len().alias("All_SNPs"))

full_snp_file = os.path.join(output_directory,f"{args.snp_name}_All_SNPs.parquet")
full_snp_df.write_parquet(full_snp_file,compression="snappy")

if args.min_table:
    
    group_min_df = pd.read_csv(args.min_table,sep="\t")
    
    assert set(group_min_df['SNP_Group'].unique()) == set(group_df['SNP_Group'].unique()), "SNP_Group columns do not match"

    min_counts =  (
        dict(zip(group_min_df['SNP_Group'], group_min_df['Min_Samples']))
    )

elif args.prop:
    min_counts = {
        group: math.ceil(float(args.prop) * len(ids))
        for group, ids in group_data.items()
    }
    
else:
    min_counts = {
        group: 1 for group, ids in group_data.items()
    }


thresh_snp_df = (
    full_snp_df
    .with_columns([
        pl.col("SNP_Group").replace_strict(min_counts).alias("Min_Count")
        ])
    .filter(pl.col("Fixed_Count") >= pl.col("Min_Count"))
).select(["contig_index","contig_position","SNP_Group","SNP_Base","Missing"])

thresh_snp_count_df = thresh_snp_df.group_by("SNP_Group").agg(pl.len().alias("Threshold_Filter"))

if args.missing:
    max_missing = args.missing
else:
    max_missing = len(sample_ids) - 1

missing_snp_df = thresh_snp_df.filter(pl.col("Missing") <= max_missing)

missing_snp_count_df = missing_snp_df.group_by("SNP_Group").agg(pl.len().alias("Missing_Filter"))

thresh_snp_file = os.path.join(output_directory,f"{args.snp_name}_Threshold_SNPs.parquet")
missing_snp_df.write_parquet(thresh_snp_file,compression="snappy")

snp_count_file = os.path.join(output_directory,f"{args.snp_name}_SNP_Counts.tsv")

(
    full_snp_count_df
    .join(thresh_snp_count_df, on="SNP_Group", how="left")
    .join(missing_snp_count_df, on="SNP_Group", how="left")
    .with_columns([
        pl.exclude("SNP_Group").fill_null(0)
    ])
    .sort("SNP_Group")
    .write_csv(snp_count_file, separator="\t")
)

if args.top:
    
    top_snp_file = os.path.join(output_directory,f"{args.snp_name}_Top_{args.top}_SNPs.parquet")

    top_snp_df = (
        thresh_snp_df
        .sort('Missing')
        .group_by("SNP_Group")
        .head(args.top)  
    )    
        
    top_snp_df.write_parquet(top_snp_file,compression="snappy")


