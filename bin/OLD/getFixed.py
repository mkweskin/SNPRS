#!/usr/bin/env python3

import polars as pl
import os
import sys
import argparse
import json
from natsort import natsorted
import pandas as pd
import math

def saveFixedSites(code_parquet, sample_list, group_name,min_count,out_path):
    
    fixed_set = {1,2,3,4,16,33,34,35,36,48}
    
    def normalize_code(col):
        return pl.when(col >= 33).then(col - 32).otherwise(col)

    if len(sample_list) == 1:
        
        sample_id = sample_list[0]
        fixed = (
            pl.scan_parquet(code_parquet)
            .select(['contig_index', 'contig_position'] + sample_id)
            .filter(pl.col(sample_id).is_in(fixed_set))
            .with_columns([pl.col(sample_id).cast(pl.Int32).alias(group_name)])
            .with_columns([normalize_code(pl.col(group_name)).alias(group_name)])
            .select(["contig_index","contig_position",group_name])
            .sort(['contig_index','contig_position'])
        ).collect(engine="streaming")
    
    else:
        
        fixed = (
            pl.scan_parquet(code_parquet)
            .select(['contig_index', 'contig_position'] + sample_list)
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
            .filter(pl.col("Fixed_Count") >= min_count)
            .with_columns([normalize_code(pl.col(group_name)).alias(group_name)])
            .select(["contig_index","contig_position",group_name])
            .sort(['contig_index','contig_position'])
            ).collect(engine="streaming")
    
    fixed.write_parquet(out_path,compression="snappy")

#####

##### Args #####

parser = argparse.ArgumentParser(description='Generate alignment from SNPRS data')

# Set path to JSON file and output directory
parser.add_argument('-j','--json_file',dest="json_file", type=str,required=True, help='Path to SNPRS joined JSON file')
parser.add_argument('-o','--out',dest="out_dir", type=str,required=True, help='Path to output directory')

# Define groups via single TSV... 
parser.add_argument('-g','--groups',dest="group_file", default = None, type=str, help='Path to TSV file with group information (Sample_ID, SNP_Group)')

# ...or a file with a list + name
parser.add_argument('-i','--id_file',dest="id_file", default = None, type=str, help='Path to text file with list of IDs')
parser.add_argument('-n','--group_name',dest="group_name", default = None, type=str, help='ID to use for SNP group if IDs provided via -i')

# Set the required proportion of samples that have data...
parser.add_argument('-p','--fixed_prop',dest="fixed_prop", default=1.0,type=float, help='Proportion of samples within a group required to call a fixed site [Default: 1.0 (100%)]')

# ... or set discrete limits via TSV
parser.add_argument('--min_table',dest="min_table", default=None,type=str, help='Path to TSV file with group count information (SNP_Group,Min_Samples)')

#####

##### Main #####

args = parser.parse_args()

output_directory = os.path.abspath(args.out_dir)

# Parse JSON
json_file = os.path.abspath(args.json_file)

with open(json_file, "r") as f:
    data = json.load(f)

sample_ids = natsorted(data["Sample_IDs"].split(","))
sample_count = len(sample_ids)

join_id = data["Join_ID"]
joined_directory = data["Joined_Directory"]
code_file = data["Code_File"]

# Parse SNP groups
if args.group_file:
    
    group_file = os.path.abspath(args.group_file)
    group_df = pd.read_csv(group_file,sep="\t")

    group_data =  (
        group_df.groupby("SNP_Group")["Sample_ID"]
        .apply(list)
        .to_dict()
    )

    group_ids = list({sid for ids in group_data.values() for sid in ids})
    assert all(item in group_ids for item in sample_ids), f"Not all samples in {args.group_file} are represented in {args.json_file}"

elif args.id_file and args.group_name:

    group_data = {}
    id_list = []
    
    with open(args.id_file, 'r') as file:
        for line in file:
            sample_id = line.strip()
            if sample_id in sample_ids:
                id_list.append(sample_id)
            else:
                sys.exit(f"{sample_id} from {args.id_file} not found in {args.json_file}")
        
    group_data[args.group_name] = id_list

else:
    sys.exit("Must specify -g or -i/-n")

# Gather or calculate the number of sites required to call a fixed site
snp_groups = list(group_data.keys())

if args.min_table:
    
    group_min_df = pd.read_csv(args.min_table,sep="\t")
    min_list = list(group_min_df['SNP_Group'].unique())
    
    assert all(item in snp_groups for item in min_list), f"Not all SNP groups in {args.min_table} are defined as groups"

    min_counts =  (
        dict(zip(group_min_df['SNP_Group'], group_min_df['Min_Samples']))
    )

else:
    min_counts = {
        group: math.ceil(float(args.fixed_prop) * len(ids))
        for group, ids in group_data.items()
    }

# Call fixed sites for each group
for i,(group, ids) in enumerate(group_data.items()):
    
    out_path = os.path.join(output_directory, f"{group}_Fixed.parquet")
    
    if os.path.exists(out_path):
        sys.exit(f"{out_path} exists...")

    group_min = min_counts[group]
        
    saveFixedSites(code_file, ids, group, group_min,out_path)






