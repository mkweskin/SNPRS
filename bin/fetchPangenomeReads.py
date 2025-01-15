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
        sys.exit("Both Forward and Reverse are empty for " + str(sample_id))
    
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
        sys.exit("Reverse read exists without forward read for " + str(sample_id))
    
    return sample_id, pangenome_group, reads

def pair_reads(read_files, subdir, forward_suffix, reverse_suffix):
    
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
        
    # Process SE reads
    single_end = [x for x in read_files if x not in left_pairs + right_pairs]

    # Create a dataframe to store read information (Sample_ID,subdir,Reads)
    read_df = pd.DataFrame(columns=['Sample_ID','Pangenome_Group','Reads'])
    
    # Add paired reads to dataframe
    for left, right in zip(left_pairs, right_pairs):
        base = str(os.path.basename(left).replace(forward_suffix,""))
        read_df = read_df.append({'Sample_ID': base, 'Pangenome_Group': str(subdir), 'Reads': ";".join([left,right])}, ignore_index=True)
    
    # Add single reads to dataframe
    for single in single_end:
        if single.endswith(forward_suffix):
            base = str(os.path.basename(single).replace(forward_suffix,""))
        else:
            base = str(os.path.basename(single).replace(read_filetype,""))
        read_df = read_df.append({'Sample_ID': base, 'Pangenome_Group': str(subdir), 'Reads': single}, ignore_index=True)

    return read_df

# Parse args
parser = argparse.ArgumentParser(description='Fetch Reads')
parser.add_argument('--read_dir', type=str, help='Path to directory containing read files or 4-column CSV with read information (Sample_ID, Pangenome_Group, Forward, Reverse)')
parser.add_argument('--read_filetype',default='fastq.gz', type=str, help='Read filetype information')
parser.add_argument('--forward_suffix',default='_1.fastq.gz', type=str, help='Forward suffix')
parser.add_argument('--reverse_suffix',default = '_2.fastq.gz', type=str, help='Reverse suffix')
args = parser.parse_args()

# Get read filetype information
read_filetype = args.read_filetype
if not read_filetype.startswith("."):
    read_filetype = "." + read_filetype
forward_suffix = args.forward_suffix
reverse_suffix = args.reverse_suffix

# Establish read_df
read_df = pd.DataFrame(columns=['Sample_ID','Pangenome_Group','Reads'])

# Check if read_dir is a csv file
read_dir = args.read_dir
if not os.path.isdir(read_dir) and not read_dir.endswith('.csv'):
    sys.exit("--read_dir is not a valid directory or .csv file: " + str(read_dir))
elif read_dir.endswith('.csv'):
    read_csv = pd.read_csv(read_dir, sep=',',header=None,dtype=str)
    if read_csv.shape[1] != 4:
        sys.exit("CSV file must have 4 columns: Sample_ID, Pangenome_Group, Forward, Reverse")
    if read_csv.iloc[0,0] == 'Sample_ID':
        read_csv = read_csv.drop(0)
    read_csv.columns = ['Sample_ID','Pangenome_Group','Forward','Reverse']
    
    for index, row in read_csv.iterrows():
        sample_id, pangenome_group, reads = parse_row(row)
        read_df = read_df.append({'Sample_ID': sample_id, 'Pangenome_Group': pangenome_group, 'Reads': reads}, ignore_index=True)
else:
    subdirs = [x for x in next(os.walk(read_dir))[1] if len(glob(os.path.join(read_dir, x, f"*{read_filetype}"))) > 0]
    unknown_reads = sorted(glob(os.path.join(read_dir, f"*{read_filetype}")))

    if len(subdirs) + len(unknown_reads) == 0:
        sys.exit("No "+read_filetype +" files detected in "+read_dir)
    
    # Process reads in subdirectories
    if len(subdirs) > 0:
        for subdir in subdirs:
            read_files = sorted(glob(os.path.join(read_dir, subdir, f"*{read_filetype}")))
            read_df = read_df.append(pair_reads(read_files, subdir, forward_suffix, reverse_suffix), ignore_index=True)
    
    # Process unknown reads
    if len(unknown_reads) > 0:
        read_df = read_df.append(pair_reads(unknown_reads, "Unknown", forward_suffix, reverse_suffix), ignore_index=True)

if read_df.empty:
    sys.exit("No reads detected from " + str(read_dir))
elif read_df['Sample_ID'].duplicated().any():
    print(read_df[read_df.duplicated(subset='Sample_ID', keep=False)])
    sys.exit("Duplicate sample IDs detected in read_df")
elif read_df['Reads'].duplicated().any():
    print(read_df[read_df.duplicated(subset='Reads', keep=False)])
    sys.exit("Duplicate reads detected in read_df")
else:
    read_df.to_csv(sys.stdout, index=False, header=False)