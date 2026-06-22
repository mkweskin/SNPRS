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


def code_to_base(x):
    base_convert_dict = { 1: "A", 2: "C", 3: "G", 4: "T" }
    return base_convert_dict.get(x, "N")

args = parse_args()

site_info = pl.scan_parquet(args.scaffold_parquet).select(['contig_index','contig_position'])
called_info = pl.scan_parquet(args.called_parquet).select(['contig_index','contig_position','base_code'])

joined = (
    site_info
    .join(called_info, on=["contig_index", "contig_position"], how="left")
    .with_columns(
        pl.col("base_code").fill_null(0)
    )
    .sort(['contig_index','contig_position'])
).collect(engine="streaming")

bases = [code_to_base(x) for x in joined["base_code"]]

sequence = "".join(bases)

sample_name = os.path.basename(args.called_parquet)
if sample_name.endswith("_Called.parquet"):
    sample_name = sample_name.replace("_Called.parquet", "")
else:
    sample_name = sample_name.replace(".parquet", "")

def wrap80(seq):
    return "\n".join(seq[i:i+80] for i in range(0, len(seq), 80))

out_path = f"{sample_name}.fasta"

with open(out_path, "w", newline="\n") as f:
    f.write(f">{sample_name}\n")
    f.write(wrap80(sequence))
    f.write("\n")