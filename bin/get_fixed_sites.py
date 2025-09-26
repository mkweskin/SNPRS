import polars as pl
import os
import glob
import pandas as pd
from typing import List
from concurrent.futures import ThreadPoolExecutor,as_completed,ProcessPoolExecutor
import pyarrow.parquet as pq
import json 
import gzip
import argparse
import sys
import psutil
from datetime import datetime
import subprocess
import time
from natsort import natsorted
import re
import shutil
import math
import numpy as np
import io
from Bio import SeqIO
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Get fixed sites from a joined SNPRS output or a subset of samples")
    
    # Data args
    parser.add_argument("--joined", dest="joined_dir", type=str, required=True,help="Path to directory containing joined output (Bases/Scaffold/Codes/Sites/Missing)")
    parser.add_argument("--group", dest="group_info", type=str, required=True,help="Path to a file containing Sample IDs to group")
    parser.add_argument("--missing", dest="missing_info", type=str, default=None,help="Path to a file containing Sample IDs allowed to have missing data")
    parser.add_argument("--name", dest="fixed_name", type=str, required=True,help="Prefix for output files")
    parser.add_argument("--fasta", dest="ref_fasta", type=str, required=True,help="Path to associated reference assembly")
    parser.add_argument("--out", dest="output_directory", type=str, required=True,help="Path to store output files and temp directory")
    
    # Filter args
    parser.add_argument("--gaps", dest="include_gaps", action="store_true",help="Include positions where 1+ samples has a gap (-) [Default: FALSE]")

    return parser.parse_args()

def get_fixed_sites(fixed_id,scaffold_file,base_file,code_file,site_file,include_gaps,group_file,missing_file,output_directory,joined_directory,ref_fasta):
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/get_fixed.py")
    
    gap_string  = "1" if include_gaps else "0"
        
    cmd = ["python", worker_script,fixed_id,scaffold_file,base_file,code_file,site_file,gap_string,group_file,missing_file,output_directory,joined_directory,ref_fasta]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to filter valid site types.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
        
##### MAIN #####

# region 00: Parse args and set up directories
args = parse_args()

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

joined_directory = os.path.abspath(args.joined_dir)
output_directory = os.path.abspath(args.output_directory)
ref_fasta = os.path.abspath(args.ref_fasta)

if not os.path.exists(joined_directory):
    sys.exit(f"${joined_directory} does not exist...")
if not os.path.exists(output_directory):
    sys.exit(f"${output_directory} does not exist...")
if not os.path.exists(ref_fasta):
    sys.exit(f"${ref_fasta} does not exist...")

# Fetch joined files
scaffold_file = list(Path(joined_directory).glob("*_Scaffold.parquet"))[0]
join_id = re.sub(r"_Scaffold$", "", scaffold_file.stem)

scaffold_file = os.path.join(joined_directory, f"{join_id}_Scaffold.parquet")
code_file = os.path.join(joined_directory, f"{join_id}_Codes.parquet")
site_file = os.path.join(joined_directory, f"{join_id}_Sites.parquet")
base_file = os.path.join(joined_directory, f"{join_id}_Bases.parquet")
called_bases_file = os.path.join(joined_directory, f"{join_id}_Called_Bases.txt")

# Filter prefix
fixed_id = str(args.fixed_name)
output_json = os.path.join(output_directory, f"{fixed_id}.json")

# Site inclusion/removal
include_gaps = args.include_gaps

# Get sample count + IDs
sample_ids = []
sample_count = 0

with open(called_bases_file, "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        sample_count += 1
        sample_id = os.path.basename(line).replace("_Called.parquet","")
        sample_ids.append(sample_id)

if sample_count < 1:
    sys.exit(f"Fewer than 1 *_Called.parquet files provided in {called_bases_file}")

# Samples to group
group_samples = []
if args.group_info:
    with open(args.group_info, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            group_samples.append(line)

# Missing data
missing_allowed = []
if args.missing_info:
    with open(args.missing_info, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            missing_allowed.append(line)

missing_not_in_group = set(missing_allowed) - set(group_samples)
assert not missing_not_in_group, f"Samples in missing_info not found in group_info: {missing_not_in_group}"

group_not_in_samples = set(group_samples) - set(sample_ids)
assert not group_not_in_samples, f"Samples in group_info not found in called bases: {group_not_in_samples}"

# Process reference FASTA
raw_records = [(rec.id, str(rec.seq)) for rec in SeqIO.parse(ref_fasta, "fasta")]
if not raw_records:
    sys.exit("No contigs found.")

raw_records = natsorted(raw_records, key=lambda x: x[0])
total_sites = sum(len(seq) for _, seq in raw_records)

contig_data = [(rec_id.strip().split()[0],len(seq)) for rec_id, seq in raw_records]
index_key = {index: (contig_id,contig_length) for index,(contig_id, contig_length) in enumerate(contig_data)}
contig_key = {contig_id: (index,contig_length) for index,(contig_id, contig_length) in enumerate(contig_data)}

# endregion

# region 01: Site Filtering
get_fixed_sites(fixed_id,scaffold_file,base_file,code_file,site_file,include_gaps,str(args.group_info),str(args.missing_info),output_directory,joined_directory,ref_fasta)

# endregion