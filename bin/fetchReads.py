#!/usr/bin/env python3

import os
import sys
from glob import glob
import argparse
import pandas as pd

# Create process to parse rows from read_csv
def parse_row(row):
    sample_id = row['Sample_ID']
    pangenome_group = row['Pangenome_Group']
    forward = row['Forward']
    reverse = row['Reverse']
    
    # Check for data
    if pd.isnull(sample_id):
        sys.exit("Sample_ID is empty for row: " + str(row))
    if pd.isnull(forward) and pd.isnull(reverse):
        sys.exit("Both Forward and Reverse are empty for row: " + str(row))
    
    # Set pangenome group to unknown if empty
    if pd.isnull(pangenome_group):
        pangenome_group = "Unknown"

    if not pd.isnull(forward):
        forward = os.path.abspath(forward)
        if not os.path.exists(forward):
            sys.exit("Forward read does not exist: " + str(forward))
    else:
        forward = []
    
    if not pd.isnull(reverse):
        reverse = os.path.abspath(reverse)
        if not os.path.exists(reverse):
            sys.exit("Reverse read does not exist: " + str(reverse))
    else:
        reverse = []
    
    if forward and reverse:
        reads = ";".join([forward,reverse])
    elif forward:
        reads = forward
    else:
        reads = reverse
    
    return sample_id, pangenome_group, reads

# Parse args
parser = argparse.ArgumentParser(description='Fetch Reads')
parser.add_argument('--read_dir', type=str, help='path to directory containing read files or 4-column CSV with read information (Sample_ID, Pangenome_Group, Forward, Reverse)')
parser.add_argument('--read_filetype',default='fastq.gz', type=str, help='read filetype information')
parser.add_argument('--forward_suffix',default='_1.fastq.gz', type=str, help='forward suffix')
parser.add_argument('--reverse_suffix',default = '_2.fastq.gz', type=str, help='reverse suffix')
args = parser.parse_args()

# Get read filetype information
read_filetype = args.read_filetype
if not read_filetype.startswith("."):
    read_filetype = "." + read_filetype
forward_suffix = args.forward_suffix
reverse_suffix = args.reverse_suffix

# Establish read_df
read_df = pd.DataFrame(columns=['Sample_ID','Pangenome_Group','Reads'])

# Check if read_dir is a tsv file
read_dir = os.path.abspath(args.read_dir)

if not os.path.isdir(read_dir) and not read_dir.endswith('.csv'):
    sys.exit("--read_dir is not a valid directory or .csv file: " + str(read_dir))
elif read_dir.endswith('.csv'):
    read_csv = pd.read_csv(read_dir, sep=',',header=None)
    if read_csv.shape[1] != 4:
        sys.exit("CSV file must have 4 columns: Sample_ID, Pangenome_Group, Forward, Reverse")
    if read_csv.iloc[0,0] == 'Sample_ID':
        read_csv = read_csv.drop(0)
    read_csv.columns = ['Sample_ID','Pangenome_Group','Forward','Reverse']
    
    # Rea
    
    
else:
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # Check if sequence files exist in directory, ignoring undetermined reads
    read_files = sorted(glob(read_dir+"/*"+read_filetype))
    read_files = [r for r in read_files if not os.path.basename(r).startswith("Undetermined_")]

    if len(read_files) == 0:
        sys.exit("No "+read_filetype +" files detected in "+read_dir)

    # Get data for paired-end and single-end reads
    left_files = [s for s in read_files if s.endswith(forward_suffix)]
    right_files = [s for s in read_files if s.endswith(reverse_suffix)]

    # Identify pairs based on file name
    left_pairs = list()
    right_pairs = list()
    paired_files = list(set([x.replace(forward_suffix, '') for x in left_files]).intersection([y.replace(reverse_suffix, '') for y in right_files]))

    for pair in paired_files:
        left_pairs.append(pair+forward_suffix)
        right_pairs.append(pair+reverse_suffix)
    single_end = [x for x in read_files if x not in left_pairs + right_pairs]

    for left in left_pairs:
        base = str(os.path.basename(left).replace(forward_suffix,"").replace(trim_name,""))
        print(",".join([base,"Paired",";".join([left,left.replace(forward_suffix,reverse_suffix)])]))

    for single in single_end:
        if single.endswith(forward_suffix):
            base = str(os.path.basename(single).replace(forward_suffix,"").replace(trim_name,""))
        else:
            base = str(os.path.basename(single).replace(read_filetype,"").replace(trim_name,""))
        print(",".join([base,"Single",single]))
