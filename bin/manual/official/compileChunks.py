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

import time
import statistics
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
from collections import defaultdict
import multiprocessing
import random
import string
import math
from datetime import datetime

def parse_args():
    parser = argparse.ArgumentParser(description="Create scaffold parquet from 1+ called base files")
    
    parser.add_argument("--parquets", dest="chunk_parquets", type=str, required=True,help="Path to directory with called base parquets, or path to a file with 1+ paths called base parquets")
    parser.add_argument("--out", dest="scaffold_file", type=str, required=True,help="Path to output file")
    parser.add_argument("--batch", dest="batch_size", type=int, default=1,help="Batch size for data processing")
    parser.add_argument("--join_id", dest="join_id", type=str, required=True,help="ID for join")
    parser.add_argument("--keep_temp", action="store_true",help="Do not delete temp files")
    return parser.parse_args()

# region 00: Parse args and set up directories
args = parse_args()
parquet_file = os.path.abspath(args.chunk_parquets)
output_parquet = os.path.abspath(args.scaffold_file)

output_directory = Path(output_parquet).parent
join_id = args.join_id

called_base_file = Path(os.path.join(output_directory,f"{join_id}_Called_Bases.txt"))
chunk_file = Path(os.path.join(output_directory,f"{join_id}_Chunk_Parquets.txt"))

with open(called_base_file, 'r') as f:
    sample_count = len(f.readlines())
    
with open(parquet_file, "r") as f:
    batch_files = [line.strip() for line in f if line.strip()]

batches = [
    batch_files[i:i + args.batch_size]
    for i in range(0, len(batch_files), args.batch_size)
]

def aggregate_batch(files):
    return (
        pl.concat([pl.scan_parquet(f) for f in files], how="vertical")
        .group_by("key")
        .sum()
        .collect(engine="streaming")
    )

batches = [
    batch_files[i:i + args.batch_size]
    for i in range(0, len(batch_files), args.batch_size)
]

acc = aggregate_batch(batches[0])

batch_count = 0
for batch in batches[1:]:
    batch_count+=args.batch_size
    df = aggregate_batch(batch)

    acc = (
        pl.concat([acc, df])
        .group_by("key")
        .sum()
    )
    
cols = ["a", "c", "g", "t", "gap"]

try:

    (
        acc
        .with_columns([
            (pl.col("key") // (2**32)).cast(pl.Int32).alias("contig_index"),
            (pl.col("key") & (2**32 - 1)).cast(pl.Int32).alias("contig_position"),
        ])
        .with_columns([
            pl.sum_horizontal([pl.col(col) >= 2 for col in cols]).alias("pi_alleles"),
            pl.sum_horizontal([pl.col(col) == 1 for col in cols]).alias("sing_count"),
            pl.sum_horizontal(pl.col(cols)).alias("fixed")
        ])
        .with_columns(
            uncovered=(pl.lit(sample_count) - pl.col("cov"))
        )
        .sort(['contig_index','contig_position'])
        .select([
            'contig_index','contig_position','cov','uncovered','fixed','het','pf','pi_alleles','a','c','g','t','gap'
        ])

        .write_parquet(output_parquet, compression="snappy")
    )
    
except:
    raise Exception("lolwut")
else:
    
    if not args.keep_temp:
            
        if chunk_file.exists():
            chunk_file.unlink()
            
        for temp_dir in output_directory.glob(f"Temp_{join_id}_*"):
            if temp_dir.is_dir():
                shutil.rmtree(temp_dir)