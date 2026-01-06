import sys
import polars as pl
import os
import json
from collections import defaultdict, Counter
from natsort import natsorted
from concurrent.futures import ProcessPoolExecutor,as_completed
import math
import argparse
import pyarrow.parquet as pq
import glob
import multiprocessing
import time
import pandas as pd

fixed_codes = {1, 2, 3, 4, 16}
het_codes = set(range(5, 16)) | set(range(17, 32))
cpu = multiprocessing.cpu_count()

def batched(iterable, batch_size):
    batch = []
    for item in iterable:
        batch.append(item)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch

def classify_batch(rows):
    counts_list = []
    codes_list = []
    for row in rows:
        counts, codes = classify_row(row)
        counts_list.append(counts)
        codes_list.append(codes)
    return counts_list, codes_list

def parse_args():
    parser = argparse.ArgumentParser(description="Create scaffold parquet from called base files")
    parser.add_argument("--chunk_dir", dest="chunk_dir", type=str, required=True,help="Path to chunk directory with scaffolded samples")
    parser.add_argument("--scaffold_parquet", dest="scaffold_parquet", type=str, required=True,help="Path to scaffold parquet")
    parser.add_argument("--chunk_tsv", dest="chunk_tsv", type=str, required=True,help="Path to TSV file with chunk information")

    return parser.parse_args()

def classify_row(row):

    uncovered_count = 0
    ploidy_fail_count = 0
    het_count = 0
    fixed_counts = {}

    for b in row:
        if b == 0:
            uncovered_count += 1
        elif b < 0:
            ploidy_fail_count += 1
        else:
            if b in het_codes:
                het_count += 1
            if b in fixed_codes:
                fixed_counts[b] = fixed_counts.get(b, 0) + 1

    nonsingleton_sample_count = 0
    nonsingleton_allele_count = 0
    singleton_count = 0

    if fixed_counts:
        nonsingleton_bases = {b: c for b, c in fixed_counts.items() if c > 1}
        singleton_count = sum(1 for c in fixed_counts.values() if c == 1)

        if nonsingleton_bases:
            nonsingleton_sample_count = sum(nonsingleton_bases.values())
            nonsingleton_alleles = frozenset(nonsingleton_bases.keys())
            row = [(32 + b) if b in nonsingleton_alleles else b for b in row]

            nonsingleton_allele_count = len(nonsingleton_alleles)

    total_missing = uncovered_count + ploidy_fail_count

    if singleton_count == 0:
        site_code = nonsingleton_allele_count
    elif nonsingleton_allele_count == 0:
        site_code = 0
    else:
        site_code = 5 + nonsingleton_allele_count

    return_count = {
        "Site_Code": site_code,
        "Missing": total_missing,
        "Uncovered": uncovered_count,
        "Filtered": ploidy_fail_count,
        "Fixed": nonsingleton_sample_count + singleton_count,
        "Hets": het_count,
        "Singletons": singleton_count,
        "Nonsingletons": nonsingleton_sample_count,
        "Nonsingleton_Alleles": nonsingleton_allele_count
    }

    return return_count, row

def wrap(row):
    return classify_row(list(row))

def compile_chunks(sorted_files, out_file):
    for f in sorted_files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Missing chunk file: {f}")

    lazy_frames = [pl.scan_parquet(f) for f in sorted_files]
    (
        pl.concat(lazy_frames, how="vertical")
        .collect(streaming=True)
        .write_parquet(out_file, compression="snappy")
    )

    for f in sorted_files:
        os.remove(f)

##### MAIN #####
args = parse_args()

chunk_dir = os.path.abspath(args.chunk_dir)
if not os.path.exists(chunk_dir):
    sys.exit(f"{chunk_dir} does not exist.")
chunk_id = os.path.basename(os.path.normpath(chunk_dir))

chunk_tsv = os.path.abspath(args.chunk_tsv)
if not os.path.exists(chunk_tsv):
    sys.exit(f"{chunk_tsv} does not exist.")

scaffold_parquet = os.path.abspath(args.scaffold_parquet)
if not os.path.exists(scaffold_parquet):
    sys.exit(f"{scaffold_parquet} does not exist.")
    
chunk_df = pd.read_csv(chunk_tsv, sep="\t")
chunk_row = chunk_df.loc[chunk_df["Chunk_ID"] == chunk_id]
if chunk_row.empty:
    raise ValueError(f"Chunk_ID {chunk_id} not found in chunk_df")

start_value = int(chunk_row["Start"].iloc[0])
stop_value = int(chunk_row["Stop"].iloc[0])

length = stop_value - start_value
if length <= 0:
    raise ValueError(f"Invalid slice length: {length} for Chunk_ID {chunk_id}")

sliced_scaffold = (
    pl.scan_parquet(scaffold_parquet)
    .slice(start_value, length)
    .collect()
)
    
output_site_file = os.path.join(chunk_dir,f"{chunk_id}_Sites.parquet")
output_code_file = os.path.join(chunk_dir,f"{chunk_id}_Codes.parquet")

if os.path.exists(output_site_file):
    sys.exit(f"{output_site_file} already exists...")
if os.path.exists(output_code_file):
    sys.exit(f"{output_code_file} already exists...")
    
parquet_files = [f for f in os.listdir(chunk_dir) if f.endswith(".parquet")]
if len(parquet_files) == 0:
    sys.exit(f"No parquet files found in {chunk_dir}...")

sample_ids = natsorted([os.path.splitext(f)[0] for f in parquet_files])
sorted_parquets = [os.path.join(chunk_dir, f"{sid}.parquet") for sid in sample_ids]

lazy_frames = [pl.scan_parquet(f) for f in sorted_parquets]
chunk_df = pl.concat(lazy_frames, how="horizontal").collect()
chunk_rows = [list(row) for row in chunk_df.iter_rows()]

num_workers = os.cpu_count()
batch_size = max(1, len(chunk_rows) // (num_workers * 4))

all_counts = []
all_codes = []

with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
    futures = [executor.submit(classify_batch, batch) for batch in batched(chunk_rows, batch_size)]
    for future in futures:
        counts_batch, codes_batch = future.result()
        all_counts.extend(counts_batch)
        all_codes.extend(codes_batch)

pl.concat([sliced_scaffold,pl.DataFrame(all_counts, orient="row")], how="horizontal").write_parquet(output_site_file, compression="snappy")
pl.concat([sliced_scaffold,pl.DataFrame(all_codes, schema=sample_ids, orient="row")], how="horizontal").write_parquet(output_code_file, compression="snappy")

print(chunk_dir)