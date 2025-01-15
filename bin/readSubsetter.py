#!/usr/bin/env python3

import os
import sys
import re
import pandas as pd
import argparse
from pathlib import Path

### Main Script ###

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genomesize', type=int, required=True, help="Genome size estimate in bp")
parser.add_argument('-c', '--coverage', type=float, default=10, help="Desired coverage")
parser.add_argument('-r', '--reads', type=str, required=True, help="Read data in CSV (Sample_ID,Pangenome_Group,Read_Location,Read_Count,Base_Count)")
parser.add_argument('-o', '--output', type=str, required=True, help="Output directory for subset reads")
args = parser.parse_args()

# Create log file
log_file = f"{args.output}/Subset_Log.txt"
with open(log_file, 'w') as log:
    log.write("SNPRS Read Subsetter Log\n")
    log.write("-------------------------------------------------------\n\n")
    log.write(f"\t- Genome size estimate: {args.genomesize} bp\n")
    log.write(f"\t- Desired coverage: {args.coverage}X\n")
    log.write(f"\t- Read Information: {os.path.abspath(args.reads)}\n")
    log.write(f"\t- Output directory: {os.path.abspath(args.output)}\n")
    log.write("\n-------------------------------------------------------\n\n")
        
# Read in read information
read_df = pd.read_csv(args.reads)
pangenome_groups = read_df['Pangenome_Group'].unique()

subset_depth = int((args.coverage * args.genomesize) / len(pangenome_groups))

with open(log_file, 'a') as log:
    log.write(f"\t- Found {len(pangenome_groups)} pangenome groups\n")
    log.write(f"\t- Requested subset depth: {subset_depth} bp per group\n")

# Prepare DataFrame
read_df['ToSample'] = 0

for group in read_df['Pangenome_Group'].unique():
    group_df = read_df[read_df['Pangenome_Group'] == group]
    total_bases = group_df['Base_Count'].sum()
    if total_bases < subset_depth:
        with open(log_file, 'a') as log:
            log.write(f"\t- NOTE: {group} samples have only {total_bases} bp. All reads will be used.\n")
        read_df.loc[read_df['Pangenome_Group'] == group, 'ToSample'] = group_df['Base_Count']
    else:
        per_dataset = subset_depth / len(group_df)
        low_datasets = group_df[group_df['Base_Count'] < per_dataset]['Sample_ID'].tolist()
        high_datasets = group_df[group_df['Base_Count'] >= per_dataset]['Sample_ID'].tolist()
        read_df.loc[read_df['Sample_ID'].isin(low_datasets), 'ToSample'] = read_df[read_df['Sample_ID'].isin(low_datasets)]['Base_Count']
        group_df = group_df[group_df['Sample_ID'].isin(high_datasets)]
        per_dataset = int((subset_depth - read_df[read_df['Sample_ID'].isin(low_datasets)]['Base_Count'].sum()) / len(group_df))
        read_df.loc[read_df['Sample_ID'].isin(high_datasets), 'ToSample'] = per_dataset

read_df.to_csv(f"{args.output}/Subset_Scheme.csv", index=False)

with open(log_file, 'a') as log:
    log.write("\n-------------------------------------------------------\n")
    log.write(f"\t- Finished!\n")
    log.write(f"\t- Subset scheme written to {args.output}/Subset_Scheme.csv\n")
    log.write("-------------------------------------------------------\n")

read_df.to_csv(sys.stdout, index=False, header=False)