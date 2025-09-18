#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path


# Parse args
parser = argparse.ArgumentParser(description='Fetch Reads')
parser.add_argument('-r','--read_count',dest="read_count", required=True,type=str, help='Path to Read_Counts.csv')
parser.add_argument('-o','--output',dest="output_directory", required=True,type=str, help='Path to directory to link reads')
parser.add_argument('--extension ',dest="read_filetype",default='fastq.gz', type=str, help='Read extension ')
parser.add_argument('--forward',dest = "forward_suffix",default='_1.fastq.gz', type=str, help='Forward suffix')
parser.add_argument('--reverse',dest = "reverse_suffix",default = '_2.fastq.gz', type=str, help='Reverse suffix')
args = parser.parse_args()

read_count_file = os.path.abspath(args.read_count)
read_count_df = pd.read_csv(read_count_file,header=None,names=['Sample_ID','Read_Count','Base_Count','Forward','Reverse'])

if not read_count_df["Sample_ID"].is_unique:
    seen = {}
    new_ids = []
    for sid in read_count_df["Sample_ID"]:
        if sid not in seen:
            seen[sid] = 0
            new_ids.append(sid)
        else:
            seen[sid] += 1
            new_ids.append(f"{sid}_{seen[sid]}")
    read_count_df["Sample_ID"] = new_ids

output_directory = os.path.abspath(args.output_directory)

for _, row in read_count_df.iterrows():
    sample_id = row["Sample_ID"]
    fwd = row["Forward"]
    rev = row["Reverse"]

    if pd.notna(fwd) and pd.notna(rev) and fwd and rev:
        target_fwd = output_directory / f"{sample_id}.{args.forward}"
        target_rev = output_directory / f"{sample_id}.{args.reverse}"
        os.symlink(Path(fwd).resolve(), target_fwd)
        os.symlink(Path(rev).resolve(), target_rev)
    elif pd.notna(fwd) and fwd:
        target_fwd = output_directory / f"{sample_id}.{args.forward}"
        os.symlink(Path(fwd).resolve(), target_fwd)
    elif pd.notna(rev) and rev:
        sys.exit("Reverse only read?")