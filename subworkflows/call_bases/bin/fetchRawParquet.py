#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import numpy as np
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
    with open(raw_parquet_data) as f:
        raw_parquet_files = [os.path.abspath(line.strip()) for line in f if line.strip() and line.strip().endswith("_Raw.parquet")]

if len(raw_parquet_files) == 0:
    sys.exit(f"No BAM files detected via {raw_parquet_data}")

raw_parquet_tuples = [(os.path.splitext(os.path.basename(raw_parquet_file))[0], raw_parquet_file) for raw_parquet_file in raw_parquet_files]

for sample, raw_parquet_file in raw_parquet_tuples:
    sample = sample[:-4]
    print(f"{sample},{raw_parquet_file}")
