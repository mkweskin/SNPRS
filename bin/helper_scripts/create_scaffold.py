import sys
import polars as pl
import os
from natsort import natsorted
import pyarrow.parquet as pq
import pandas as pd
import argparse
import shutil
import subprocess
import csv

def parse_args():
    parser = argparse.ArgumentParser(description="Create scaffold parquet from called base files")
    
    parser.add_argument("--called_bases", dest="called_base_file", type=str, required=True,help="File with paths to 2+ Called_Bases parquet")
    parser.add_argument("--join_id", dest="join_id", type=str, required=True,help="Prefix for output files")
    parser.add_argument("--out_dir", dest="output_directory", type=str, default=None,help="Path to output Parquet files [Default: cwd]")
    parser.add_argument("--mem_factor", dest="memory_factor", type=int, default=1,help="Chunk size = CPU*mem_factor")
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

def save_scaffold_parquet(sorted_parquet_paths, output_parquet):

    scaffold = (
        pl.scan_parquet(sorted_parquet_paths[0])
        .select(['contig_index','contig_position'])
        .collect()
    )
    
    if len(sorted_parquet_paths) > 1:

        for parquet in sorted_parquet_paths[1:]:

            next_df = (
                pl.scan_parquet(parquet)
                .select(['contig_index','contig_position'])
                .collect()
            )

            scaffold = pl.concat([scaffold, next_df], how="vertical").unique()

    scaffold.sort(["contig_index", "contig_position"]).write_parquet(output_parquet,compression="snappy")
       
# region 00: Parse args and set up directories
args = parse_args()

join_id = str(args.join_id)
mem_factor = int(args.memory_factor)

# Output directory
if args.output_directory is None:
    output_directory = os.getcwd()
else:
    output_directory = os.path.abspath(args.output_directory)
if not os.path.exists(output_directory):
    sys.exit(f"{output_directory} (--out_dir) does not exist")

# Temp directory
temp_directory = os.path.join(output_directory,f"Temp_{join_id}")
if not os.path.exists(temp_directory):
    sys.exit(f"{temp_directory} does not exist")

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
    save_scaffold_parquet(sorted_parquet_paths,output_parquet)
except Exception as e:
    sys.exit(f"ERROR creating scaffold: {e}")

# endregion

# region 03: Save chunk TSV
chunk_tsv_file = os.path.join(temp_directory,"Chunk_Info.tsv")

row_count = pq.ParquetFile(output_parquet).metadata.num_rows

max_chunks = os.cpu_count() * (mem_factor)
n_chunks = min(max_chunks, row_count)
chunk_size = max(1, (row_count + n_chunks - 1) // n_chunks)

chunk_info = []

for i in range(n_chunks):

    chunk_id = f"Chunk_{i}"
    chunk_dir = os.path.join(temp_directory,chunk_id)

    start = i * chunk_size
    stop = min(start + chunk_size, row_count)

    if start >= stop:
        continue
    
    os.mkdir(chunk_dir)
    chunk_info.append([f"Chunk_{i}",chunk_dir,start,stop])

df = pd.DataFrame(chunk_info, columns=["Chunk_ID","Chunk_Directory", "Start", "Stop"])
df.to_csv(chunk_tsv_file, sep="\t", index=False)

print(",".join([output_parquet,chunk_tsv_file]))

# endregion