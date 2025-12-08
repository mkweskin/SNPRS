import sys
import polars as pl
import os
from natsort import natsorted
import pyarrow.parquet as pq
import argparse
import shutil
import subprocess
import csv

def parse_args():
    parser = argparse.ArgumentParser(description="Create scaffold parquet from called base files")
    
    parser.add_argument("--called_bases", dest="called_base_file", type=str, required=True,help="File with paths to 2+ Called_Bases parquet")
    parser.add_argument("--join_id", dest="join_id", type=str, required=True,help="Prefix for output files")
    parser.add_argument("--out_dir", dest="output_directory", type=str, default=None,help="Path to output Parquet files [Default: cwd]")
    return parser.parse_args()

def fetch_base_parquets(file_path):

    with open(file_path, "r") as f:
        paths = [line.strip() for line in f if line.strip()]

    for path in paths:
        if not os.path.exists(path):
            print(f"Error: Parquet file '{path}' not found.")
            sys.exit(1)
    
    if len(paths) < 1:
        print("Error: You must provide at least one called base parquet files.")
        sys.exit(1)

    return [os.path.abspath(path) for path in paths]

def get_scaffold_sites(temp_directory,sorted_parquet_paths):
    
    valid_sites = [0, 1, 3, 4, 6]
    
    tmp_file = os.path.join(temp_directory, "Temp.parquet")
    stage_file = os.path.join(temp_directory, "Stage.parquet")

    (
        pl.scan_parquet(sorted_parquet_paths[0])
        .filter(pl.col("type").is_in(valid_sites))
        .select(["contig_index", "contig_position"])
        .unique()
        .sort(['contig_index','contig_position'])
        .sink_parquet(tmp_file,compression = "snappy")
    )

    if len(sorted_parquet_paths) > 1:
        
        for path in sorted_parquet_paths[1:]:
            
            lazy_scaffold = pl.scan_parquet(tmp_file)
            
            lazy_new = (
                pl.scan_parquet(path)
                .filter(pl.col("type").is_in(valid_sites))
                .select(["contig_index", "contig_position"])
            )

            (
                pl.concat([lazy_scaffold, lazy_new])
                .unique()
                .sink_parquet(stage_file)
            )

            shutil.move(stage_file, tmp_file)

    return tmp_file

def save_scaffold_parquet(output_parquet,temp_parquet):
    
    (
        pl.scan_parquet(temp_parquet)
        .sort(['contig_index','contig_position'])
        .sink_parquet(output_parquet,compression = "snappy")
    )

# region 00: Parse args and set up directories
args = parse_args()

join_id = str(args.join_id)

# Output directory
if args.output_directory is None:
    output_directory = os.getcwd()
else:
    output_directory = os.path.abspath(args.output_directory)
if not os.path.exists(output_directory):
    sys.exit(f"{output_directory} (--out_dir) does not exist")

# Temp directory
temp_directory = os.path.join(output_directory,f"Temp_{join_id}")
if os.path.exists(temp_directory):
    shutil.rmtree(temp_directory)
os.mkdir(temp_directory)

# Called Bases
called_base_file = os.path.abspath(args.called_base_file)
if not os.path.exists(called_base_file):
    sys.exit(f"{called_base_file} (--called_bases) does not exist")

# Output parquet
output_parquet = os.path.join(output_directory,f"{join_id}_Scaffold.parquet")
if os.path.exists(output_parquet):
    sys.exit(f"{output_parquet} already exists")

# endregion

# region 01: Fetch called base files

called_base_files = fetch_base_parquets(called_base_file)
sample_id_paths = [(os.path.basename(path).replace("_Called.parquet", ""), path) for path in called_base_files]
sorted_path_sample_pairs = natsorted(sample_id_paths, key=lambda x: x[0])
sorted_samples, sorted_parquet_paths = zip(*sorted_path_sample_pairs)

# endregion

# region 02: Save scaffold file

try:
    temp_parquet = get_scaffold_sites(temp_directory,sorted_parquet_paths)
    
    if len(sorted_parquet_paths) > 1:
        save_scaffold_parquet(output_parquet, temp_parquet)
    else:
        shutil.move(temp_parquet, output_parquet)

finally:
    shutil.rmtree(temp_directory)

# endregion