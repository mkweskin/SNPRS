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
    
    parser.add_argument("--called", dest="called_bases", type=str, required=True,help="Path to directory with called base parquets, or path to a file with 1+ paths called base parquets")
    parser.add_argument("--output", dest="output_directory", type=str, required=True,help="Output directory")
    parser.add_argument("--name", dest="output_id", type=str, required=True,help="Output ID")
    parser.add_argument("--batch", dest="batch_size", type=int, default=1,help="Batch size for data processing")
    
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

def generate_random_string(length):
    characters = string.ascii_letters + string.digits
    return ''.join(random.choices(characters, k=length))

# region 00: Parse args and set up directories
args = parse_args()

# Batch size (Default 1)
batch_size = int(args.batch_size)

# Output directory
output_directory = os.path.abspath(args.output_directory)
if not os.path.exists(output_directory):
    sys.exit(f"{output_directory} does not exist...")

# Output ID
output_id = str(args.output_id)

# Output parquet
output_parquet = os.path.join(output_directory,f"{output_id}_Scaffold.parquet")
if os.path.exists(output_parquet):
    sys.exit(f"{output_parquet} already exists")
    
# Temp directory
now = datetime.now()
timestamp = now.strftime("%Y_%m_%d_%H_%M_%S")
temp_directory = os.path.join(output_directory,f"Temp_{output_id}_{timestamp}")
if os.path.exists(temp_directory):
    sys.exit(f"{temp_directory} exists...")
else:
    os.mkdir(temp_directory)
    sys.stderr.write(f"\n\t- Created TEMP directory at {temp_directory} ...\n")

# Called Bases
called_bases = os.path.abspath(args.called_bases)

if os.path.isdir(called_bases):
    called_base_files = glob.glob(os.path.join(called_bases, "*_Called.parquet"))
elif os.path.isfile(called_bases):
    called_base_files = fetch_base_parquets(called_bases)
else:
    raise ValueError(f"No valid called parquets found via --called ({called_bases})")

called_count = len(called_base_files)
num_batches = math.ceil(called_count / batch_size)

sys.stderr.write(f"\t- Read in and located {called_count} called base files from {called_bases} ...\n")

batch_count = 0
batch_files = []
for i in range(0, called_count, batch_size):
    batch = called_base_files[i:i + batch_size]

    rand = generate_random_string(10)
    batch_path = os.path.join(temp_directory, f"{rand}.txt")

    with open(batch_path, "w") as f:
        f.write("\n".join(batch)+"\n")
    
    batch_count+=1
    batch_files.append(batch_path)
    

sys.stderr.write(f"\t- Broke data up into {batch_count} batches for processing (batch size: {args.batch_size})\n")

print("\n".join(batch_files))