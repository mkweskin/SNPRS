#!/usr/bin/env python3
import sys
import polars as pl
import numpy as np
import json
import struct
import os
import argparse


fixed_codes = {1, 2, 3, 4, 16}    

a_codes   = {1,5,21,6,22,7,23,11,27,12,28,13,29,15,17,31}
c_codes   = {2,5,21,8,24,9,25,11,27,12,28,14,30,15,18,31}
g_codes   = {3,6,22,8,24,10,26,11,27,13,29,14,30,15,19,31}
t_codes   = {4,7,23,9,25,10,26,12,28,13,29,14,30,15,20,31}
gap_codes = {16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}
   
het_codes = list(set(range(5, 16)) | set(range(17, 32)))

def parse_args():
    parser = argparse.ArgumentParser(description="Score chunk of scaffold")    
    parser.add_argument("--txt", dest="txt_file", type=str, required=True,help="Path to text file with paths to called base files")
    return parser.parse_args()

args = parse_args()
chunk_txt = os.path.abspath(args.txt_file)
root, ext = os.path.splitext(chunk_txt)
output_parquet = root + ".parquet"
with open(chunk_txt, "r") as f:
    batch_files = [line.strip() for line in f if line.strip()]
    
(
    pl.concat(
        [
            pl.scan_parquet(p)
            .select(["contig_index", "contig_position", "base_code"])
            .with_columns([
                (pl.col("contig_index").cast(pl.UInt64) * (2**32) +
                    pl.col("contig_position").cast(pl.UInt64))
                .alias("key")])
            .select("key", "base_code")
            for p in batch_files
        ],
        how="vertical"
    )
    .with_columns([
        (pl.col("base_code").is_in(het_codes)).cast(pl.UInt32).alias("het"),
        (pl.col("base_code") < 0).cast(pl.UInt32).alias("pf"),
        (pl.col("base_code") == 1).cast(pl.UInt32).alias("a"),
        (pl.col("base_code") == 2).cast(pl.UInt32).alias("c"),
        (pl.col("base_code") == 3).cast(pl.UInt32).alias("g"),
        (pl.col("base_code") == 4).cast(pl.UInt32).alias("t"),
        (pl.col("base_code") == 16).cast(pl.UInt32).alias("gap"),
    ])
    .group_by("key")
    .agg([
        pl.len().alias("cov"),
        pl.sum("het").alias("het"),
        pl.sum("pf").alias("pf"),
        pl.sum("a").alias("a"),
        pl.sum("c").alias("c"),
        pl.sum("g").alias("g"),
        pl.sum("t").alias("t"),
        pl.sum("gap").alias("gap"),
    ])
    .sink_parquet(output_parquet,compression="snappy")
)

print(output_parquet)