#!/usr/bin/env python3

import os
import sys
import argparse
from pathlib import Path

# Parse args

parser = argparse.ArgumentParser(description='Fetch Raw Parquets')
parser.add_argument('-p','--parquet_data',dest="parquet_data", type=str, help='Path to directory containing parquet files or a file with a list of parquet files')
args = parser.parse_args()

raw_parquet_data = os.path.abspath(args.parquet_data)

raw_parquet_files = []
if os.path.isdir(raw_parquet_data):
    raw_parquet_files = [str(f.resolve()) for f in Path(raw_parquet_data).glob("*_Raw.parquet")]
elif os.path.isfile(raw_parquet_data):

    if raw_parquet_data.endswith("_Raw.parquet"):
        raw_parquet_files = [os.path.abspath(raw_parquet_data)]

    else:
        with open(raw_parquet_data, "r", encoding="utf-8") as f:
            raw_parquet_files = [
                os.path.abspath(line.strip())
                for line in f
                if line.strip().endswith("_Raw.parquet")
            ]

else:
    raise FileNotFoundError(f"{raw_parquet_data} does not exist")
if len(raw_parquet_files) == 0:
    sys.exit(f"No *_Raw.parquet files detected via {raw_parquet_data}")

raw_parquet_tuples = [(os.path.splitext(os.path.basename(raw_parquet_file))[0], raw_parquet_file) for raw_parquet_file in raw_parquet_files]

for sample, raw_parquet_file in raw_parquet_tuples:
    sample = sample[:-4]
    print(f"{sample},{raw_parquet_file}")
