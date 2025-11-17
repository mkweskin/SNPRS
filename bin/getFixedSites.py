#!/usr/bin/env python3
import pandas as pd
import os
import sys
from natsort import natsorted
from collections import Counter
import argparse
import polars as pl
import pyarrow.parquet as pq
import shutil
import json

def parse_args():
    parser = argparse.ArgumentParser(description="Get fixed sites from 1+ _Called.parquet files")
    
    # Data args
    parser.add_argument("--out", dest="output_directory", type=str, required=True,help="Path to save output")
    parser.add_argument("--fixed_id", dest="fixed_id", type=str, required=True,help="Prefix for output files")
    parser.add_argument("--called_bases", dest="called_bases", type=str, required=True,help="Path to file with 1+ called base paths")
    parser.add_argument("--no_gaps", dest="no_gaps", action = "store_true", help="Do not include deletions in fixed sites")
    parser.add_argument("--missing", dest="missing", type=float, default=-1,help="If -1, 1 sample required [Default]. If 0 or >= 1, max number of samples allowed with missing data. If between 0 - 1, the minimum proportion of samples with data required.")
    return parser.parse_args()

def fetch_base_parquets(file_path):

    with open(file_path, "r") as f:
        paths = [line.strip() for line in f if line.strip()]

    for path in paths:
        if not os.path.exists(path):
            print(f"Error: Parquet file '{path}' not found.")
            sys.exit(1)
        
    if len(paths) == 0:
        sys.exit("No paths provided by --called_bases")

    return [os.path.abspath(path) for path in paths]

def save_singleton(output_directory, fixed_id, no_gaps,sorted_samples, sorted_parquet_paths):

    output_parquet = f"{output_directory}/{fixed_id}.parquet"
    output_json = f"{output_directory}/{fixed_id}.json"

    if no_gaps:
        valid_sites = [0]
        gap_string = "No"
    else:
        valid_sites = [0, 1]
        gap_string = "Yes"

    (
        pl.scan_parquet(sorted_parquet_paths[0])
        .filter(pl.col("type").is_in(valid_sites))
        .select([
            "contig_index",
            "contig_position",
            pl.col("final_base").alias(fixed_id),
        ])
        .unique()
        .sort(["contig_index", "contig_position"])
        .with_columns(pl.lit(1).alias("Count"))
        .sink_parquet(output_parquet)
    )

    fixed_count = pq.ParquetFile(output_parquet).metadata.num_rows

    if fixed_count == 0:
        sys.exit(f"No fixed sites detected for {sorted_samples[0]}")

    fixed_info = {
        "Samples": sorted_samples[0],
        "Called Base Files": sorted_parquet_paths[0],
        "Missing Allowed": "0",
        "Include Gaps": gap_string,
        "Fixed Site Count": fixed_count,
    }

    with open(output_json, "w", encoding="utf-8") as f:
        json.dump(fixed_info, f, indent=4)


def get_scaffold_sites(temp_directory,valid_sites,sorted_parquet_paths):
    
    tmp_file = os.path.join(temp_directory, "Temp.parquet")
    stage_file = os.path.join(temp_directory, "Stage.parquet")

    (
        pl.scan_parquet(sorted_parquet_paths[0])
        .filter(pl.col("type").is_in(valid_sites))
        .select(["contig_index", "contig_position"])
        .unique()
        .sink_parquet(tmp_file)
    )

    for path in sorted_parquet_paths[1:]:
        lazy_scaffold = pl.scan_parquet(tmp_file)
        lazy_new = (
            pl.scan_parquet(path)
            .filter(pl.col("type").is_in(valid_sites))
            .select(["contig_index", "contig_position"])
        )

        (
            pl.concat([lazy_scaffold, lazy_new])
            .unique()
            .sink_parquet(stage_file)
        )

        shutil.move(stage_file, tmp_file)

    return tmp_file

def save_site_parquet(temp_directory, fixed_id, temp_parquet):
    
    output_parquet = f"{temp_directory}/{fixed_id}_Sites.parquet"

    (
        pl.scan_parquet(temp_parquet)
        .sort(['contig_index','contig_position'])
        .sink_parquet(output_parquet,compression = "snappy")
    )

    return output_parquet

def save_fixed_parquet(output_directory,fixed_id,site_parquet,sorted_samples,sorted_parquet_paths,valid_sites,max_missing):

    gap_string = "Yes" if 1 in valid_sites else "No"

    final_parquet = f"{output_directory}/{fixed_id}.parquet"
    output_json = f"{output_directory}/{fixed_id}.json"

    lazy_joined = pl.scan_parquet(site_parquet)

    min_allowed = len(sorted_samples) - int(max_missing)

    for sample_id,path in zip(sorted_samples,sorted_parquet_paths):
        
        lazy_sample = (
            pl.scan_parquet(path)
            .filter(pl.col("type").is_in(valid_sites))
            .select(['contig_index','contig_position',pl.col('final_base').alias(sample_id)])
        )
        
        lazy_joined = (
            lazy_joined
            .join(lazy_sample,on=['contig_index','contig_position'],how="left")
        )

    
    fixed_rows = (
        lazy_joined 
        .fill_null("?")
        .with_columns([
            pl.concat_list([pl.when(pl.col(c) != "?").then(pl.col(c)).otherwise(None) for c in sorted_samples])
            .alias("nonzero_values")
        ])
        .with_columns([
            pl.col("nonzero_values").list.drop_nulls().list.unique().alias("unique_values")
        ])
        .filter(pl.col("unique_values").list.len() == 1)

        .with_columns([
            pl.col("unique_values").list.first().alias(fixed_id),
            pl.col("nonzero_values").list.drop_nulls().list.len().alias("Count")
        ])
        .filter(pl.col("Count") >= min_allowed)
        .select(['contig_index','contig_position',fixed_id,"Count"])
        .unique()
        .sort(['contig_index','contig_position'])
        .sink_parquet(final_parquet)

    )

    sample_string = ';'.join(sorted_samples)
    called_base_string = ';'.join(sorted_parquet_paths)
    fixed_count = pq.ParquetFile(final_parquet).metadata.num_rows
    
    fixed_info = {
        "Samples": sample_string,
        "Called Base Files": called_base_string,
        "Missing Allowed": f"{max_missing}",
        "Include Gaps": gap_string,
        "Fixed Site Count": fixed_count,
    }

    with open(output_json, "w", encoding="utf-8") as f:
        json.dump(fixed_info, f, indent=4)

##### MAIN #####

# region 00: Parse args

args = parse_args()

fixed_id = args.fixed_id

output_directory = os.path.abspath(args.output_directory)
if not os.path.exists(output_directory):
    sys.exit(f"{output_directory} (--out) does not exist")
temp_directory = os.path.abspath(f"{output_directory}/Temp")
os.mkdir(temp_directory)

called_base_file = os.path.abspath(args.called_bases)
if not os.path.exists(called_base_file):
    sys.exit(f"{called_base_file} (--called_bases) does not exist")

called_base_files = fetch_base_parquets(called_base_file)

# Sample Data
sample_id_paths = [(os.path.basename(path).replace("_Called.parquet", ""), path) for path in called_base_files]
sorted_path_sample_pairs = natsorted(sample_id_paths, key=lambda x: x[0])
sorted_samples, sorted_parquet_paths = zip(*sorted_path_sample_pairs)
sample_count = len(sorted_samples)


if sample_count == 1:
    save_singleton(output_directory,fixed_id,args.no_gaps,sorted_samples,sorted_parquet_paths)

else:

    # Missing data
    if float(args.missing) == -1:
        max_missing = sample_count - 1
    elif float(args.missing) == 0:
        max_missing = 0
    elif 0 < float(args.missing) < 1:
        min_present = int(float(args.missing) * sample_count)
        max_missing = sample_count - min_present
    else:
        max_missing = int(args.missing)

    # Gaps
    if args.no_gaps:
        valid_sites = [0]
        gap_string = "Include Gaps: No"
    else:
        valid_sites = [0,1]
        gap_string = "Include Gaps: Yes"

    try:
        temp_parquet = get_scaffold_sites(temp_directory,valid_sites,sorted_parquet_paths)
        site_parquet = save_site_parquet(temp_directory, fixed_id, temp_parquet)
        fixed_parquet = save_fixed_parquet(output_directory,fixed_id,site_parquet,sorted_samples,sorted_parquet_paths,valid_sites,max_missing)

    finally:
        shutil.rmtree(temp_directory)

