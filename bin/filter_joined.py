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
    parser = argparse.ArgumentParser(description="Filter SNPRS output")
    
    # Data args
    parser.add_argument("--joined", dest="joined_dir", type=str, required=True,help="Path to directory containing joined output (Bases/Scaffold/Codes/Sites/Missing)")
    parser.add_argument("--fasta", dest="ref_fasta", type=str, required=True,help="Path to associated reference assembly")
    parser.add_argument("--out", dest="output_directory", type=str, required=True,help="Path to store output files and temp directory")
    parser.add_argument("--name", dest="filter_name", type=str, required=True,help="Output prefix [Default: <ANALYSIS_NAME>_{timestamp}]")
    
    # Filter args
    parser.add_argument("--types", dest="site_types", type=str, default='btqp',help="String of single letter codes for sites requested: F/f: Fixed; B/b: Biallelic; T/t: Triallelic; Q/q: Quadallelic; P/p: Pentallellic; S/s: Singleton-only")
    parser.add_argument("--gaps", dest="include_gaps", action="store_true",help="Include positions where 1+ samples has a gap (-) [Default: FALSE]")
    parser.add_argument("--het", dest="include_hets", action="store_true",help="Include positions where 1+ samples has a heterozygous base call [Default: FALSE]")
    parser.add_argument("--invalid", dest="include_invalid", action="store_true",help="Include positions where 1+ samples has an invalid base call [Default: FALSE]")
    parser.add_argument("--nosing", dest="no_singletons", action="store_true",help="Do not include sites if any sample has a singleton allele [Default: FALSE]")
    parser.add_argument("--missing", dest="missing", type=float, default=None,help="If >= 1, max number of samples allowed with missing data. If < 1, the minimum proportion of samples with data required. [Default: Estimate from data]")

    return parser.parse_args()

def gross_exclusion(site_file,include_gaps,include_hets,remove_singletons,include_invalid,temp_directory):
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/gross_exclusion.py")
    
    gap_string  = "1" if include_gaps else "0"
    het_string  = "1" if include_hets else "0"
    invalid_string = "1" if include_invalid else "0"
    sing_string = "0" if remove_singletons else "1"
        
    cmd = ["python", worker_script,site_file,gap_string,het_string,invalid_string,sing_string,temp_directory]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to filter valid site types.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    output = result.stdout.decode("utf-8").strip().split(",")
    return output

def filter_site_types(site_file,pass_gross_exclusion_file,valid_site_types,temp_directory):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/filter_valid_types.py")
    
    fixed_string  = "1" if "f" in valid_site_types else "0"
    bi_string  = "1" if "b" in valid_site_types else "0"
    tri_string  = "1" if "t" in valid_site_types else "0"
    quad_string  = "1" if "q" in valid_site_types else "0"
    pent_string = "1" if "p" in valid_site_types else "0"
    sing_string = "1" if "s" in valid_site_types else "0"
    
    cmd = ["python", worker_script,site_file,pass_gross_exclusion_file,fixed_string,bi_string,tri_string,quad_string,pent_string,sing_string,temp_directory]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to filter valid site types.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    output = result.stdout.decode("utf-8").strip().split(",")
    return output
 
def get_max_missing(site_file,pass_site_type_file,temp_directory):
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/estimate_missing.py")
    
    cmd = ["python", worker_script,site_file,pass_site_type_file,temp_directory]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to filter valid site types.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    output = result.stdout.decode("utf-8").strip()
    return int(output)

def filter_missing(site_file,pass_site_type_file,max_missing,temp_directory):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/filter_missing.py")
    

    cmd = ["python", worker_script,site_file,pass_site_type_file,str(max_missing),temp_directory]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to filter valid site types.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    output = result.stdout.decode("utf-8").strip().split(",")
    return output

def save_filtered_parquet(raw_parquet,new_parquet,row_file):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/save_filtered_join.py")
    

    cmd = ["python", worker_script,raw_parquet,new_parquet,row_file]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to filter valid site types.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

def make_alignment(base_file,row_file,output_fasta):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/sites2align.py")
    

    cmd = ["python", worker_script,base_file,row_file,output_fasta]
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

temp_directory = os.path.join(output_directory, f"SNPRS_Temp_{timestamp}")
os.mkdir(temp_directory)

# Fetch joined files
scaffold_file = list(Path(joined_directory).glob("*_Scaffold.parquet"))[0]
join_id = re.sub(r"_Scaffold$", "", scaffold_file.stem)

scaffold_file = os.path.join(joined_directory, f"{join_id}_Scaffold.parquet")
code_file = os.path.join(joined_directory, f"{join_id}_Codes.parquet")
site_file = os.path.join(joined_directory, f"{join_id}_Sites.parquet")
base_file = os.path.join(joined_directory, f"{join_id}_Bases.parquet")
missing_file = os.path.join(joined_directory, f"{join_id}_Missing.tsv")
called_bases_file = os.path.join(joined_directory, f"{join_id}_Called_Bases.txt")

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

if sample_count < 2:
    sys.exit(f"Fewer than 2 *_Called.parquet files provided in {called_bases_file}")

# Filter prefix
filter_id = str(args.filter_name)
output_fasta = os.path.join(output_directory,f"{filter_id}_aln.fasta")
new_scaffold_file = os.path.join(output_directory, f"{filter_id}_Scaffold.parquet")
new_code_file = os.path.join(output_directory, f"{filter_id}_Codes.parquet")
new_site_file = os.path.join(output_directory, f"{filter_id}_Sites.parquet")
new_base_file = os.path.join(output_directory, f"{filter_id}_Bases.parquet")
new_missing_file = os.path.join(output_directory, f"{filter_id}_Missing.tsv")
new_called_bases_file = os.path.join(output_directory, f"{filter_id}_Called_Bases.txt")

# Site Types
requested_site_types = sorted(set(re.sub(r'[^A-Za-z]', '', args.site_types).lower()))

if len(requested_site_types) < 1:
    sys.exit("Error: no valid site types provided.")

invalid = set(requested_site_types) - set("fbtqps")

if invalid:
    sys.exit(f"Error: invalid site type(s) [Allowed: fbtqps]: {', '.join(sorted(invalid))}")
    
valid_site_types = "".join(requested_site_types)
if args.no_singletons and "s" in valid_site_types:
    sys.exit("Error: --nosing specified, but singleton sites were requested via 's' in --site_types")

# Site inclusion/removal
include_gaps = args.include_gaps
include_hets = args.include_hets
include_invalid = args.include_invalid
remove_singletons = args.no_singletons

# Missing data
if not args.missing:
    missing_mode = "estimate"
    max_missing = np.nan
elif int(args.missing) == 0:
    max_missing = 0
elif float(args.missing) < 1:
    missing_mode = "prop"
    min_present = int(float(args.missing) * sample_count)
    max_missing = sample_count - min_present
elif int(args.missing) >= 1:
    missing_mode = "value"
    max_missing = int(args.missing)

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

pass_gross_exclusion_file,pass_gross_exclusion_count = gross_exclusion(site_file,include_gaps,include_hets,remove_singletons,include_invalid,temp_directory)
pass_site_type_file,pass_site_type_count = filter_site_types(site_file,pass_gross_exclusion_file,valid_site_types,temp_directory)

if missing_mode == "estimate":
    max_missing = get_max_missing(site_file,pass_site_type_file,temp_directory)

pass_missing_file,pass_missing_count = filter_missing(site_file,pass_site_type_file,max_missing,temp_directory)

# endregion

# region 02: Save filtered output 
save_filtered_parquet(scaffold_file,new_scaffold_file,pass_missing_file)
save_filtered_parquet(code_file,new_code_file,pass_missing_file)
save_filtered_parquet(site_file,new_site_file,pass_missing_file)
save_filtered_parquet(base_file,new_base_file,pass_missing_file)

# Links for now
os.symlink(missing_file, new_missing_file)
os.symlink(called_bases_file, new_called_bases_file)

# endregion

# region 03: Save Alignment
make_alignment(base_file,pass_missing_file,output_fasta)

# endregion

# region 04: Save Summary
filtering_info = {
    "Joined_Directory":joined_directory,
    "Ref_FASTA": ref_fasta,
    "Sample_IDs":sample_ids,
    "Sample_Count":str(sample_count),
    "Output_Directory":output_directory,
    "Site_Types":valid_site_types,
    "Gaps_Included":str(include_gaps),
    "Hets_Included":str(include_hets),
    "Invalid_Included":str(include_invalid),
    "Singletons_Removed":str(remove_singletons),
    "Missing_Arg":str(args.missing),
    "Max_Missing":str(max_missing),
    "Final_Site_Count":str(pass_missing_count)
}

output_json = os.path.join(output_directory, f"{filter_id}.json")
with open(output_json, "w", encoding="utf-8") as f:
    json.dump(filtering_info, f, indent=4)

# endregion