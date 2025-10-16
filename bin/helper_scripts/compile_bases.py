import sys
import polars as pl
import os
from natsort import natsorted
import pyarrow.parquet as pq
import argparse
import shutil
import subprocess
import csv
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="Combine all Scaffolded_* files into a single base parquet")
    
    parser.add_argument("--join_id", dest="join_id", type=str, required=True,help="Prefix for output files")
    parser.add_argument("--out_dir", dest="output_directory", type=str, required=True,help="Path to output directory [Default: cwd]")

    return parser.parse_args()

def compile_base_files(sorted_files,output_parquet):

    row_counts = []
    for f in sorted_files:
        row_count = pq.ParquetFile(f).metadata.num_rows
        row_counts.append(row_count)

    if len(set(row_counts)) != 1:
        raise ValueError(f"Not all files have the same number of rows: {row_counts}")
    
    lazy_frames = [pl.scan_parquet(f) for f in sorted_files]
    
    (
        pl.concat(lazy_frames, how="horizontal")
        .collect(streaming=True)
        .write_parquet(output_parquet, compression="snappy")
    )
    
    for f in sorted_files:
        os.remove(f)

args = parse_args()
join_id = args.join_id

output_directory = os.path.abspath(args.output_directory)
if not os.path.exists(output_directory):
    sys.exit(f"{output_directory} (--out_dir) does not exist...")

output_parquet = os.path.join(output_directory,f"{join_id}_Bases.parquet")
if os.path.exists(output_parquet):
    sys.exit(f"{output_parquet} already exists...")

input_parquet_pattern = os.path.join(output_directory, "Scaffolded_*.parquet")
input_parquets = glob.glob(input_parquet_pattern)

if len(input_parquets) == 0:
    sys.exit(f"No files matching 'Scaffolded_*.parquet' found in {output_directory}...")

sorted_files = natsorted(input_parquets,key=lambda f: os.path.basename(f).replace("Scaffolded_", "").replace(".parquet", ""))

compile_base_files(sorted_files,output_parquet)