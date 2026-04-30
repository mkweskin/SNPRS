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
import glob

fixed_codes = {1, 2, 3, 4, 16}
het_codes = set(range(5, 16)) | set(range(17, 32))

def parse_args():
    parser = argparse.ArgumentParser(description="Create scaffold parquet from 1+ called base files")
    
    parser.add_argument("--called", dest="called_bases", type=str, required=True,help="Path to directory with called base parquets, or path to a file with 1+ paths called base parquets")
    parser.add_argument("--out", dest="scaffold_file", type=str, required=True,help="Output file name")
    parser.add_argument("--batch", dest="batch_size", type=int, default=1000,help="Batch size for data processing")
    
    return parser.parse_args()

def fetch_base_parquets(file_path):

    with open(file_path, "r") as f:
        paths = [os.path.abspath(line.strip())
                 for line in f
                 if line.strip()]

    if not paths:
        sys.exit("Error: You must provide at least one called base parquet file.")

    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        sys.exit("Error: Missing parquet files:\n" + "\n".join("  " + m for m in missing))

    return paths

def save_scaffold_parquet(parquet_paths, output_parquet, batch_size):

    seen = set()

    for i in range(0, len(parquet_paths), batch_size):
        
        batch = parquet_paths[i:i + batch_size]

        df = pl.concat(
            [pl.scan_parquet(p).select(['contig_index', 'contig_position'])
             for p in batch],
            how="vertical"
        ).unique().collect(engine="streaming")

        for row in df.rows():
            key = (row[0], row[1])
            seen.add(key)

    scaffold = pl.DataFrame(
        list(seen),
        schema=[
        ("contig_index", pl.Int32),
        ("contig_position", pl.Int32)],
        orient="row"
    ).sort(["contig_index", "contig_position"])

    scaffold.write_parquet(output_parquet, compression="snappy")
    

# region 00: Parse args and set up directories
args = parse_args()

# Output parquet
output_parquet = os.path.abspath(args.scaffold_file)
if os.path.exists(output_parquet):
    sys.exit(f"{output_parquet} already exists")

# Called Bases
called_bases = os.path.abspath(args.called_bases)

if os.path.isdir(called_bases):
    called_base_files = glob.glob(os.path.join(called_bases, "*_Called.parquet"))
elif os.path.isfile(called_bases):
    called_base_files = fetch_base_parquets(called_bases)
else:
    raise ValueError(f"No valid called parquets found via --called ({called_bases})")
try:
    save_scaffold_parquet(called_base_files,output_parquet,args.batch_size)
except Exception:
    raise

# endregion