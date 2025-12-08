#!/usr/bin/env python3

import os
import sys
import argparse
from pathlib import Path

# Parse args

parser = argparse.ArgumentParser(description='Fetch Called Bases')
parser.add_argument('-p','--parquet_data',dest="parquet_data", type=str, help='Path to directory containing parquet files or a file with a list of parquet files')
args = parser.parse_args()

called_parquet_data = os.path.abspath(args.parquet_data)

called_parquet_files = []
if os.path.isdir(called_parquet_data):
    called_parquet_files = [str(f.resolve()) for f in Path(called_parquet_data).glob("*_Called.parquet")]
elif os.path.isfile(called_parquet_data):
    with open(called_parquet_data) as f:
        called_parquet_files = [os.path.abspath(line.strip()) for line in f if line.strip() and line.strip().endswith("_Called.parquet")]

if len(called_parquet_files) == 0:
    sys.exit(f"No parquet files detected via {called_parquet_data}")

called_parquet_tuples = [(os.path.splitext(os.path.basename(called_parquet_file))[0], called_parquet_file) for called_parquet_file in called_parquet_files]

for sample, called_parquet_file in called_parquet_tuples:
    sample = sample[:-7]
    print(f"{sample},{called_parquet_file}")
