#!/usr/bin/env python3

import polars as pl
import os
import sys
import argparse
import json
from natsort import natsorted
import pandas as pd
import math

def snpSubtractor(code_file, focal_id, fixed_parquet, subtract_samples):
    
    lazy_focal = pl.scan_parquet(fixed_parquet)

    lazy_subtract = pl.scan_parquet(code_file).select(['contig_index', 'contig_position'] + subtract_samples)
            
    degenerate_map = {
        1:  {1,5,21,6,22,7,23,11,27,12,28,13,29,15,17,31},
        2:  {2,5,21,8,24,9,25,11,27,12,28,14,30,15,18,31},
        3:  {3,6,22,8,24,10,26,11,27,13,29,14,30,15,19,31},
        4:  {4,7,23,9,25,10,26,12,28,13,29,14,30,15,20,31},
        16: {16,17,18,29,20,21,22,23,24,25,26,27,28,29,30,31}
    }
        
    bg_sets = {}

    for allele, codes in degenerate_map.items():
        bg_sets[allele] = (
            lazy_subtract
            .filter(
                pl.any_horizontal([
                    pl.col(s).is_in(list(codes))
                    for s in subtract_samples
                ])
            )
            .select(["contig_index", "contig_position"])
        )

    filtered_snps = []

    for allele in [1, 2, 3, 4, 16]:
        bg_subset = bg_sets[allele]

        allele_filtered = (
            lazy_focal
            .filter(pl.col(focal_id) == allele)
            .join(bg_subset,
                  on=["contig_index", "contig_position"],
                  how="anti")
        )

        filtered_snps.append(allele_filtered)

    result = pl.concat(filtered_snps).collect(engine="streaming")
    return result

#####


##### Args #####

parser = argparse.ArgumentParser(description='Generate alignment from SNPRS data')
parser.add_argument('-j','--json_file',dest="json_file", type=str,required=True, help='Path to SNPRS joined JSON file')
parser.add_argument('-g','--groups',dest="group_file", required = True, type=str, help='Path to TSV file with group information (Sample_ID, SNP_Group)')

parser.add_argument('-o','--out',dest="out_dir", type=str,required=True, help='Path to output directory')
parser.add_argument('-n','--name',dest="snp_name", required=True,type=str, help='Prefix for output files')

parser.add_argument('-c','--covered_prop',dest="covered_prop", default=1.0,type=float, help='Proportion of total samples with positive base call required to call a covered site [Default: 1.0 (100%)]')
parser.add_argument('-s','--snp_prop',dest="snp_prop", default=1.0,type=float, help='Proportion of samples within a group required to call a SNP [Default: 1.0 (100%)]')
parser.add_argument('--min_table',dest="min_table", default=None,type=str, help='Path to TSV file with group count information (SNP_Group,Min_Samples)')

#####

##### Main #####

args = parser.parse_args()

output_directory = os.path.abspath(args.out_dir)

# Parse JSON
json_file = os.path.abspath(args.json_file)

with open(json_file, "r") as f:
    data = json.load(f)

join_id = data["Join_ID"]
joined_directory = data["Joined_Directory"]
sample_ids = natsorted(data["Sample_IDs"].split(","))
sample_count = len(sample_ids)
scaffold_file = data["Scaffold_File"]
code_file = data["Code_File"]
site_file = data["Site_File"]
sample_summary_file = data["Sample_Summary_File"]
site_count_file = data["Site_Count_File"]

lazy_codes = pl.scan_parquet(code_file)
lazy_sites = pl.scan_parquet(site_file)
missing_info = lazy_sites.select(['contig_index','contig_position','Missing']).collect(engine="streaming")

# Parse SNP groups
group_file = os.path.abspath(args.group_file)
group_df = pd.read_csv(group_file,sep="\t")

group_data =  (
    group_df.groupby("SNP_Group")["Sample_ID"]
    .apply(list)
    .to_dict()
)

all_ids = list({id for ids in group_data.values() for id in ids})

# Get sites that have base calls in --covered_prop of samples
assert args.covered_prop <= 1.0, "--covered_prop must be 1.0 or less"
cov_required = math.ceil(float(args.covered_prop) * sample_count)
max_missing = sample_count - cov_required

cov_sites = (
    lazy_sites
    .select(['contig_index','contig_position','Missing'])
    .filter(pl.col("Missing") <= max_missing)
    .with_columns([
        pl.lit("SNPRS_Covered").alias("SNP_Group"),
        pl.lit(0).alias("SNP_Base"),
        (sample_count - pl.col("Missing")).alias("Fixed_Count")
    ])
    .select(['contig_index','contig_position','SNP_Group','SNP_Base','Fixed_Count'])
).collect(engine="streaming")

# Get sites with coverage in sufficient samples to call a SNP
if args.min_table:
    
    group_min_df = pd.read_csv(args.min_table,sep="\t")
    
    assert set(group_min_df['SNP_Group'].unique()) == set(group_df['SNP_Group'].unique()), "SNP_Group columns do not match"

    min_counts =  (
        dict(zip(group_min_df['SNP_Group'], group_min_df['Min_Samples']))
    )

else:
    min_counts = {
        group: math.ceil(float(args.snp_prop) * len(ids))
        for group, ids in group_data.items()
    }

snp_results = []

for i,(group, ids) in enumerate(group_data.items()):
    
    print(f"Starting {group}...\n")
    
    non_focal = list(set(all_ids) - set(ids))
    min_count = min_counts[group]
    
    sp_snp_df = snpSubtractor(code_file,group,ids,min_count,non_focal)

    snp_results.append(sp_snp_df.with_columns(pl.lit(group).alias("SNP_Group"),
    pl.col(group).alias("SNP_Base")).select(['contig_index','contig_position','SNP_Group','SNP_Base','Fixed_Count']))





full_snp_df = pl.concat(snp_results).join(missing_info,on=["contig_index","contig_position"],how="left")

full_snp_count_df = full_snp_df.group_by("SNP_Group").agg(pl.len().alias("All_SNPs"))

full_snp_file = os.path.join(output_directory,f"{args.snp_name}_All_SNPs.parquet")
full_snp_df.write_parquet(full_snp_file,compression="snappy")


    
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



