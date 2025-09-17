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
from collections import defaultdict, Counter


def parse_args():
    parser = argparse.ArgumentParser(description="Process pileup file into parquet")
    
    parser.add_argument("-b", dest="called_bases", type=str, required=True,help="File with paths to 2+ Called_Bases parquet")
    parser.add_argument("-o", dest="output_directory", type=str, default=None,help="Path to output Parquet files [Default: cwd")
    parser.add_argument("-n", dest="name", type=str, default=None,help="Prefix for output files [Default: SNPPRS_{timestamp}]")
    
    return parser.parse_args()

def extract_sample_metadata(parquet_path):
    
    metadata = pq.ParquetFile(parquet_path).schema_arrow.metadata
    if metadata is None:
        raise ValueError(f"{parquet_path} has no metadata.")
    meta = {k.decode(): v.decode() for k, v in metadata.items()}
    
    row = {
        "Sample_ID": meta.get("sample_id"),
        "Pileup_File": meta.get("pileup_file"),
        "Sample_Parquet_File": meta.get("Sample_Parquet_File"),
        "Called_Base_File":parquet_path,
        "Ref_JSON": meta.get("ref_json"),
        "Percent_Covered": float(meta.get("percent_covered", 0.0)),
    }
    # Median of coverage_quartiles
    coverage_q = json.loads(meta.get("coverage_quartiles", "{}"))
    row["Median_Read_Coverage"] = float(coverage_q.get("50%", 0.0))

    qc_params = json.loads(meta.get("QC_Parameters", "{}"))
    row.update(qc_params)

    called_base_counts = json.loads(meta.get("Called_Base_Counts", "{}"))
    for k, v in called_base_counts.items():
        row[f"Called_{k}"] = v

    return row

def print_memory_usage():
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 * 1024)  # in MiB
    print(f"🔍 Memory usage: {mem:.2f} MiB")
    
def initiate_scaffold_chunks(parquet_path_file,output_directory,temp_directory,analysis_name):
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/create_scaffold.py")
    cmd = ["python",worker_script, parquet_path_file, output_directory,temp_directory,analysis_name]
    result = subprocess.run(cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True)
    
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to compile parquet.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
    try:
        output = json.loads(result.stdout)
    except json.JSONDecodeError:
        raise

    return output['scaffold_parquet'],output["chunk_indexes"],output["sorted_samples"],output['sorted_parquets']

def chunk_bases(scaffold_parquet,sample_parquet,sample_id,chunk_json,temp_directory):
    
    temp_sample_dir = os.path.join(temp_directory,sample_id)
    os.mkdir(temp_sample_dir)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/chunk_bases.py")
    
    cmd = ["python",worker_script,scaffold_parquet,sample_parquet,sample_id,chunk_json,temp_sample_dir]
    result = subprocess.run(cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True)
    
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to compile parquet.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
    try:
        output = json.loads(result.stdout)
    except json.JSONDecodeError:
        raise
    return output['chunk_files']
    
def populate_chunks(chunk_file,temp_directory):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/score_sites.py")
    
    cmd = ["python", worker_script,chunk_file,temp_directory]
    result = subprocess.run(cmd,
        stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True
        )
    
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to compile parquet.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
    try:
        output = json.loads(result.stdout)
    except json.JSONDecodeError:
        raise

    return output["site_file"],output['code_file'],output['missing_data'],output['chunk_id']

def compile_bases(base_parquet,base_file):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/compile_bases.py")
    cmd = ["python", worker_script, base_parquet,base_file]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to compile parquet.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

def compile_codes(output_file, file_list):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/compile_codes.py")
    cmd = ["python", worker_script, output_file, file_list]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to compile parquet.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

def compile_sites(output_file, file_list):
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/compile_sites.py")
    cmd = ["python", worker_script, output_file, file_list]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to compile parquet.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
                
##### MAIN #####
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# region 00: Parse args and set up directories
args = parse_args()

analysis_name = args.name
if analysis_name is None:
    analysis_name = f"SNPRS_{timestamp}"

parquet_path_file = args.called_bases

if args.output_directory is None:
    output_directory = os.getcwd()
else:
    output_directory = os.path.abspath(args.output_directory)

temp_directory = os.path.join(output_directory, f"SNPRS_Temp_{timestamp}")
parent_directory = os.path.dirname(output_directory)

if not os.path.exists(output_directory):
    if not os.path.exists(parent_directory):
        sys.exit(f"{parent_directory} does not exist, cannot create output folder at {output_directory}")
    os.mkdir(output_directory)
if not os.path.exists(temp_directory):
    os.mkdir(temp_directory)

# endregion

# region 01: Create scaffold file and chunks
start_time = datetime.now()
print(f"[{start_time}] Starting scaffolding...")

scaffold_parquet,chunk_json,sorted_samples,sorted_parquets = initiate_scaffold_chunks(parquet_path_file,output_directory,temp_directory,analysis_name)

end_time = datetime.now()
elapsed = end_time - start_time
print(f"[{end_time}] Scaffolding completed in {elapsed}\n")
# endregion

# region 02: Chunk bases for each sample
start_time = datetime.now()

chunk_files = []
for sample_parquet, sample_id in zip(sorted_parquets, sorted_samples):
    start_chunk_time = datetime.now()
    print(f"[{start_time}] Starting chunking for {sample_id}...")
    chunk_files.append(
        chunk_bases(scaffold_parquet, sample_parquet, sample_id, chunk_json, temp_directory)
    )
    end_chunk_time = datetime.now()
    chunk_elapsed = end_chunk_time - start_chunk_time
    print(f"[{end_chunk_time}] Chunking completed in {chunk_elapsed}\n")

organized_chunks = [list(items) for items in zip(*chunk_files)]

chunk_files = []
for i, chunk_list in enumerate(organized_chunks):
    chunk_file = os.path.join(temp_directory, f"Base_Chunks_{i}.json")
    chunk_files.append(chunk_file)
    with open(chunk_file, "w") as f:
        json.dump(chunk_list, f, indent=2)
        
end_time = datetime.now()
elapsed = end_time - start_time
print(f"[{end_time}] All chunking completed in {elapsed}\n")
# endregion

# region 03: Get codes and site info for each chunk
start_time = datetime.now()
print(f"[{start_time}] Starting scoring...")

site_list = []
code_list = []
missing_list = []
with ThreadPoolExecutor() as executor:
    # Submit all jobs at once
    futures = {
        executor.submit(populate_chunks, chunk_file, temp_directory): chunk_file for chunk_file in chunk_files
    }

    # Wait for results as they finish
    for future in as_completed(futures):
        i = futures[future]
        try:
            result = future.result()
            site_list.append(result[0])
            code_list.append(result[1])
            missing_list.append(result[2])
            print(f"Finished scoring chunk {result[3]}")
        except Exception as e:
            print(f"Error in chunk {i}: {e}")

end_time = datetime.now()
elapsed = end_time - start_time
print(f"[{end_time}] All scoring completed in {elapsed}\n")

# endregion

# region 04: Process missing data
summary = defaultdict(Counter)
for chunk in missing_list:
    for m, codes in chunk.items():
        summary[int(m)].update({int(k): int(v) for k, v in codes.items()})

# Define code labels
code_labels = {
    0: "Other",
    1: "Pure_Fixed",
    2: "Pure_Biallelic",
    3: "Pure_Triallelic",
    4: "Pure_Quadallelic",
    5: "Pure_Pentallelic",
    6: "Fixed_wSingleton",
    7: "Biallelic_wSingleton",
    8: "Triallelic_wSingleton",
    9: "Quadallelic_wSingleton",
    10: "Pentallelic_wSingleton"
}

# Save missing data
rows = [
    {"Missing": m, **{code_labels.get(c, "Other"): v for c, v in counts.items()}}
    for m, counts in summary.items()
]

df = pl.DataFrame(rows).fill_null(0)
nonzero_cols = ["Missing"] + [
    col for col in df.columns if col != "Missing" and df[col].sum() > 0
]

missing_file = os.path.join(output_directory,f"{analysis_name}_Missing.tsv")
df.select(nonzero_cols).sort("Missing").to_pandas().to_csv(missing_file, index=False,sep='\t')

# endregion

# region 05: Save final datasets
base_file = os.path.join(temp_directory,"All_Bases.txt")
site_file = os.path.join(temp_directory,"All_Sites.txt")
code_file = os.path.join(temp_directory,"All_Codes.txt")

with open(base_file, "w") as f:
    f.write("\n".join(chunk_files))
with open(site_file, "w") as f:
    f.write("\n".join(site_list))
with open(code_file, "w") as f:
    f.write("\n".join(code_list))

base_parquet = os.path.join(output_directory,f"{analysis_name}_Bases.parquet")
site_parquet = os.path.join(output_directory,f"{analysis_name}_Sites.parquet")
code_parquet = os.path.join(output_directory,f"{analysis_name}_Codes.parquet")

try:
    # Save full base parquet
    start_time = datetime.now()
    print(f"[{start_time}] Starting bases...")
    
    compile_bases(base_parquet,base_file)

    end_time = datetime.now()
    elapsed = end_time - start_time
    print(f"[{end_time}] Compiled bases in {elapsed}\n")
    
    # Save full site parquet
    start_time = datetime.now()
    print(f"[{start_time}] Starting sites...")

    compile_sites(site_parquet,site_file)
    end_time = datetime.now()
    elapsed = end_time - start_time
    print(f"[{end_time}] Compiled sites in {elapsed}\n")
    
    # Save full code parquet
    start_time = datetime.now()
    print(f"[{start_time}] Starting codes...")

    compile_codes(code_parquet,code_file)
    end_time = datetime.now()
    elapsed = end_time - start_time
    print(f"[{end_time}] Compiled codes in {elapsed}\n")
    
    shutil.rmtree(temp_directory)

except Exception as e:
    print(f"Error in final save: {e}")

# endregion    

