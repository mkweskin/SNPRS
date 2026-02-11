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
    parser = argparse.ArgumentParser(description="Based off a scaffold parquet, get base information for a sample")
    
    parser.add_argument("--called_bases", dest="called_base_file", type=str, required=True,help="Path to a called base parquet")
    parser.add_argument("--join_id", dest="join_id", type=str, required=True,help="Prefix for output files")
    parser.add_argument("--scaffold_parquet", dest="scaffold_parquet", type=str, required=True,help="Path to scaffold parquet")
    parser.add_argument("--chunk_tsv", dest="chunk_tsv", type=str, required=True,help="Path to TSV file with chunk information")
    parser.add_argument("--out_dir", dest="output_directory", type=str, required=True,help="Path to output directory [Default: cwd]")

    return parser.parse_args()

# region 00: Parse args
args = parse_args()

join_id = str(args.join_id)

# Called base file
called_base_file = os.path.abspath(args.called_base_file)
if not os.path.exists(called_base_file):
    sys.exit(f"{called_base_file} (--called_bases) does not exist")

schema = pq.read_schema(called_base_file)
metadata_bytes = schema.metadata or {}
og_metadata = {k.decode("utf-8"): v.decode("utf-8") for k, v in metadata_bytes.items()}
sample_id = og_metadata['sample_id']

# Output directory
output_directory = os.path.abspath(args.output_directory)
if not os.path.exists(output_directory):
    sys.exit(f"{output_directory} (--out_dir) does not exist")

# Temp directory
temp_directory = os.path.join(output_directory,f"Temp_{join_id}")
if not os.path.exists(temp_directory):
    sys.exit(f"{temp_directory} does not exist")

# Scaffold parquet file
scaffold_parquet = os.path.abspath(args.scaffold_parquet)
if not os.path.exists(scaffold_parquet):
    sys.exit(f"{scaffold_parquet} (--scaffold_parquet) does not exist")
    
# Chunk TSV file
chunk_tsv = os.path.abspath(args.chunk_tsv)
if not os.path.exists(chunk_tsv):
    sys.exit(f"{chunk_tsv} (--chunk_tsv) does not exist")

# endregion

# region 01: Scaffold sample and save chunks
lazy_scaffold = pl.scan_parquet(scaffold_parquet)

lazy_called = (
    pl.scan_parquet(called_base_file)
    .select(['contig_index','contig_position',pl.col("base_code").alias(sample_id)])
    )

sample_column = (
    lazy_scaffold
    .join(lazy_called,on=["contig_index","contig_position"],how="left")
    .with_columns(pl.col(sample_id).fill_null(0))
    .select(sample_id)
    .cast({sample_id: pl.Int8 })
).collect()

chunk_df = pd.read_csv(chunk_tsv, sep="\t")

for _, row in chunk_df.iterrows():

    chunk_id = row["Chunk_ID"]
    chunk_directory = row["Chunk_Directory"]
    start = int(row["Start"])
    stop = int(row["Stop"])
    chunk_file = os.path.join(chunk_directory, f"{sample_id}.parquet")
    
    if os.path.exists(chunk_file):
        os.remove(chunk_file)
        
    try:
        (
            sample_column
            .slice(start, stop - start)
            .write_parquet(chunk_file, compression="snappy")
        )
    
    except Exception as e:
        sys.exit(f"ERROR chunking sample {sample_id} for chunk {chunk_id}: {e}")
    
print(chunk_tsv)

# endregion
