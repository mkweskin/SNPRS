#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path

def pair_reads(read_dir, read_filetype, forward_suffix, reverse_suffix):
    
    read_dir = Path(read_dir)
    files = [f for f in read_dir.iterdir() if f.is_file() and f.name.endswith(read_filetype)]
    
    pe_samples = {}
    se_samples = {}
    used_files = []
    
    for f in files:
        stem = str(f)
        if stem in used_files:
            continue
        
        if stem.endswith(forward_suffix):
            sample_base = stem[:-len(forward_suffix)]
            sample_id = os.path.basename(sample_base)
            forward = f"{sample_base}{forward_suffix}" 
            reverse = f"{sample_base}{reverse_suffix}" 
            if os.path.exists(reverse):
                pe_samples[sample_id] = (os.path.abspath(forward),os.path.abspath(reverse))
                used_files.append(forward)
                used_files.append(reverse)
            else:
                se_samples[sample_id] = (os.path.abspath(forward),None)
                used_files.append(forward)

        elif stem.endswith(reverse_suffix):
            sample_base = stem[:-len(reverse_suffix)]
            sample_id = os.path.basename(sample_base)
            forward = f"{sample_base}{reverse_suffix}"
            reverse = f"{sample_base}{reverse_suffix}"
            se_samples[sample_id] = (os.path.abspath(reverse),None)
        else:
            sample_base = stem[:-len(read_filetype)]
            sample_id = os.path.basename(sample_base)
            se_samples[sample_id] = (os.path.abspath(stem),None)
    
    paired_reads = []
    for sample_id, (fwd, rev) in pe_samples.items():
        paired_reads.append((sample_id, fwd, rev))
    
    for sample_id, (fwd, rev) in se_samples.items():
        paired_reads.append((sample_id, fwd, rev))
    return pd.DataFrame(paired_reads,columns = ['Sample_ID','Forward','Reverse'])

# Parse args

parser = argparse.ArgumentParser(description='Fetch Reads')
parser.add_argument('-d','--dir',dest="read_dir", type=str, help='Path to directory containing read files or a file with a list of directories')
parser.add_argument('-e','--extension ',dest="read_filetype",default='fastq.gz', type=str, help='Read extension ')
parser.add_argument('-f','--forward',dest = "forward_suffix",default='_1.fastq.gz', type=str, help='Forward suffix')
parser.add_argument('-r','--reverse',dest = "reverse_suffix",default = '_2.fastq.gz', type=str, help='Reverse suffix')
args = parser.parse_args()

# Get read filetype information
read_filetype = args.read_filetype
if not read_filetype.startswith("."):
    read_filetype = "." + read_filetype

forward_suffix = args.forward_suffix
reverse_suffix = args.reverse_suffix

read_dir = os.path.abspath(args.read_dir)
if os.path.isdir(read_dir):
    read_dirs = [read_dir]

elif os.path.isfile(read_dir):
    with open(read_dir) as f:
        read_dirs = [os.path.abspath(line.strip()) for line in f if line.strip()]
else:
    sys.exit(f"Error: {read_dir} is neither a directory nor a file")

read_data = []
for dir in read_dirs:
    read_data.append(pair_reads(dir,read_filetype,forward_suffix,reverse_suffix))

read_df = pd.concat(read_data).reset_index(drop=True)

if not read_df["Sample_ID"].is_unique:
    seen = {}
    new_ids = []
    for sid in read_df["Sample_ID"]:
        if sid not in seen:
            seen[sid] = 0
            new_ids.append(sid)
        else:
            seen[sid] += 1
            new_ids.append(f"{sid}_{seen[sid]}")
    read_df["Sample_ID"] = new_ids

read_df.to_csv(sys.stdout, sep=",", index=False, header=False)
