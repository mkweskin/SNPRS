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

    if called_parquet_data.endswith("_Called.parquet"):
        called_parquet_files = [os.path.abspath(called_parquet_data)]

    else:
        with open(called_parquet_data, "r", encoding="utf-8") as f:
            called_parquet_files = [
                os.path.abspath(line.strip())
                for line in f
                if line.strip().endswith("_Called.parquet")
            ]

else:
    raise FileNotFoundError(f"{called_parquet_data} does not exist")

if len(called_parquet_files) == 0:
    sys.exit(f"No parquet files detected via {called_parquet_data}")

called_parquet_tuples = [(os.path.splitext(os.path.basename(called_parquet_file))[0], called_parquet_file) for called_parquet_file in called_parquet_files]

for sample, called_parquet_file in called_parquet_tuples:
    sample = sample[:-7]
    print(f"{sample},{called_parquet_file}")
