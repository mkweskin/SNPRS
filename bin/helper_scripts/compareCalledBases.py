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
import itertools

def parse_args():
    parser = argparse.ArgumentParser(description="Get pairwise comparisons between called base files")
    
    # Data args
    parser.add_argument("--out", dest="output_directory", type=str, required=True,help="Path to save output")
    parser.add_argument("--compare_id", dest="compare_id", type=str, required=True,help="Prefix for output files")
    parser.add_argument("--called_bases", dest="called_bases", type=str, required=True,help="Path to file with 1+ called base paths")
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

def compareCalledBases(sample_pair, file_pair):
    s1, s2 = sample_pair
    f1, f2 = file_pair

    called_1 = (
        pl.scan_parquet(f1)
        .filter(pl.col('type').is_in([0,1]))
        .select(['contig_index','contig_position',pl.col('final_base').alias(s1)])
    )

    het_count_1 =  (
        pl.scan_parquet(f1)
        .filter(pl.col('type').is_in([3,4]))
        .select(pl.len()).collect().item()
    )

    called_2 = (
        pl.scan_parquet(f2)
        .filter(pl.col('type').is_in([0,1]))
        .select(['contig_index','contig_position',pl.col('final_base').alias(s2)])
    )

    het_count_2 =  (
        pl.scan_parquet(f1)
        .filter(pl.col('type').is_in([3,4]))
        .select(pl.len()).collect().item()
    )

    sample_1_count = called_1.select(pl.len()).collect().item()
    sample_2_count = called_2.select(pl.len()).collect().item()

    sample_1_het_ratio = 0 if sample_1_count == 0 else het_count_1/(sample_1_count + het_count_1)
    sample_2_het_ratio = 0 if sample_2_count == 0 else het_count_2/(sample_2_count + het_count_2)

    stats = (
        called_1
        .join(called_2,on=['contig_index', 'contig_position'],how="full")
        .with_columns([
            ((pl.col(s1).is_not_null()) & (pl.col(s2).is_not_null())).alias("both"),
            ((pl.col(s1).is_not_null()) & (pl.col(s2).is_not_null()) & (pl.col(s1) == pl.col(s2))).alias("match"),
        ])
        .select([
            pl.col("both").cast(int).sum().alias("Cocalled_Count"),
            pl.col("match").cast(int).sum().alias("Match_Count"),
        ])
        .collect()
        .to_dict(as_series=False)
    )

    summary = pl.DataFrame({
        "Sample_1": [s1],
        "Sample_2": [s2],
        "Sample_1_Percent_Het": [sample_1_het_ratio],
        "Sample_2_Percent_Het": [sample_2_het_ratio],
        "Cocalled_Count": [stats["Cocalled_Count"][0]],
        "Match_Count": [stats["Match_Count"][0]]
    }).with_columns(
        pl.when(pl.col("Cocalled_Count") > 0)
          .then(pl.col("Match_Count") / pl.col("Cocalled_Count"))
          .otherwise(0.0)
          .alias("Match_Ratio")
    )

    return summary

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

compare_id = args.compare_id

output_directory = os.path.abspath(args.output_directory)
if not os.path.exists(output_directory):
    sys.exit(f"{output_directory} (--out) does not exist")

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
    sys.exit("compareCalledBases.py requires 2+ called base files as input")

else:
    
    sample_pairs = list(itertools.combinations(sorted_samples, 2))
    file_pairs = list(itertools.combinations(sorted_parquet_paths, 2))
    comparison_list = []

    paired = list(zip(sorted_samples, sorted_parquet_paths))

    comparison_list = []
    for (s1, f1), (s2, f2) in itertools.combinations(paired, 2):
        comparison_list.append(compareCalledBases((s1, s2), (f1, f2)))

    final_comparisons = pl.concat(comparison_list, how="vertical") if comparison_list else pl.DataFrame()

    long_df = (
        pl.concat([
            final_comparisons
            .select([
                pl.col("Sample_1").alias("Sample"),
                pl.col("Sample_1_Percent_Het").alias("Percent_Het"),
                pl.col("Cocalled_Count"),
                pl.col("Match_Ratio")]),
            final_comparisons.select([
                pl.col("Sample_2").alias("Sample"),
                pl.col("Sample_2_Percent_Het").alias("Percent_Het"),
                pl.col("Cocalled_Count"),
                pl.col("Match_Ratio")])])
    )

    sample_summary = (
        long_df
        .group_by("Sample")
        .agg([
            pl.first("Percent_Het").alias("Percent_Het"),
            pl.median("Cocalled_Count").alias("Median_Cocalled"),
            pl.median("Match_Ratio").alias("Median_Match_Ratio")
        ])
        .sort("Sample")
    )

    print(sample_summary)