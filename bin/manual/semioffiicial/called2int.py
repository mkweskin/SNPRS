#!/usr/bin/env python3
import sys
import polars as pl
import numpy as np
import json
import struct
import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Score chunk of scaffold")    
    parser.add_argument("--called", dest="called_parquet", type=str, required=True,help="Path to called base parquet")
    parser.add_argument("--scaffold", dest="scaffold_parquet", type=str, required=True,help="Path to scaffold base parquet")
    return parser.parse_args()

args = parse_args()

sample_name = os.path.basename(args.called_parquet)
if sample_name.endswith("_Called.parquet"):
    sample_name = sample_name.replace("_Called.parquet", "")
else:
    sample_name = sample_name.replace(".parquet", "")

out_path = f"{sample_name}.bin"

    
site_info = pl.scan_parquet(args.scaffold_parquet).select(['contig_index','contig_position'])
called_info = pl.scan_parquet(args.called_parquet).select(['contig_index','contig_position','base_code'])

joined = (
    site_info
    .join(called_info, on=["contig_index", "contig_position"], how="left")
    .with_columns(
        pl.when(pl.col("base_code").is_null())
        .then(0)
        .when(pl.col("base_code").is_in([1,2,3,4,16]))
        .then(pl.col("base_code"))
        .otherwise(-1)
        .cast(pl.Int8)
        .alias("base_code")
    )
    .sort(["contig_index", "contig_position"])
    .select("base_code")
)

sink = open(out_path, "wb")

for batch in joined.collect(engine="streaming").iter_slices(n_rows=1_000_000):
    arr = batch.to_numpy().astype(np.int8, copy=False)
    arr.tofile(sink)

sink.close()