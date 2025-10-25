import polars as pl
import os
import glob
import pandas as pd
from typing import List
from concurrent.futures import ThreadPoolExecutor,as_completed,ProcessPoolExecutor
import pyarrow.parquet as pq
import json 
import gzip
import argparse
import sys
import psutil
from datetime import datetime
import subprocess
import time
from natsort import natsorted
import re
import shutil
import math
import numpy as np
import io
from Bio import SeqIO
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Filter SNPRS output")
    
    # Data args
    parser.add_argument("--join_dir", dest="joined_dir", type=str, required=True,help="Path to directory containing joined output (Bases/Scaffold/Codes/Sites/Missing)")
    parser.add_argument("--out_dir", dest="output_directory", type=str, required=True,help="Path to store output files and temp directory")
    parser.add_argument("--filter_id", dest="filter_id", type=str, required=True,help="Output prefix [Default: <ANALYSIS_NAME>_{timestamp}]")
    parser.add_argument("--join_id", dest="join_id", type=str, required=True,help="Output prefix [Default: <ANALYSIS_NAME>_{timestamp}]")
    parser.add_argument("--fasta", dest="ref_fasta", type=str, required=True,help="Path to associated reference assembly")
    
    # Filter args
    parser.add_argument("--types", dest="site_types", type=str, default='btqp',help="String of single letter codes for sites requested: F/f: Fixed; B/b: Biallelic; T/t: Triallelic; Q/q: Quadallelic; P/p: Pentallellic; S/s: Singleton-only; U/u: Unique singletons (Fixed + Singletons)")
    parser.add_argument("--no_gaps", dest="remove_gaps", action="store_true",help="Remove positions where any sample has a gap (-) [Default: FALSE]")
    parser.add_argument("--no_sing", dest="no_singletons", action="store_true",help="Remove positions where any sample has a singleton [Default: FALSE; Overridden if S/S or U/u included in --types]")

    parser.add_argument("--het", dest="include_hets", action="store_true",help="Include positions where any sample has a heterozygous base call [Default: FALSE]")
    parser.add_argument("--invalid", dest="include_invalid", action="store_true",help="Include positions where any sample has an invalid base call [Default: FALSE]")
    parser.add_argument("--missing", dest="missing", type=float, default=None,help="If 0 or >= 1, max number of samples allowed with missing data. If between 0 -1, the minimum proportion of samples with data required. [Default: Estimate from data]")

    return parser.parse_args()

def gross_exclusion(site_file,remove_gaps,include_hets,remove_singletons,include_invalid,temp_directory):
    
    gross_exclusion_file = os.path.join(temp_directory, "Pass_Gross_Exclusion.parquet")

    sites = pl.scan_parquet(site_file)
    row_count = pq.ParquetFile(site_file).metadata.num_rows
    row_numbers = pl.LazyFrame({"row_nr": list(range(row_count))})

    expr = pl.lit(True)

    if not include_hets:
        expr &= pl.col("Hets") == 0
    if not include_invalid:
        expr &= pl.col("Filtered") == 0
    if remove_singletons:
        expr &= pl.col("Singletons") == 0
    if remove_gaps:
        expr &= pl.col("Het_Gap") == 0
        expr &= pl.col("Singleton_Gap") == 0
        expr &= pl.col("Nonsingleton_Gap") == 0

    pl.concat([sites,row_numbers],how='horizontal').filter(expr).select("row_nr").collect(streaming=True).write_parquet(gross_exclusion_file)
    
    return gross_exclusion_file

def filter_site_types(site_file,valid_site_types,temp_directory):
    
    valid_sites_file = os.path.join(temp_directory, "Pass_Valid_Sites.parquet")
    sites = pl.scan_parquet(site_file)
    row_count = pq.ParquetFile(site_file).metadata.num_rows
    row_numbers = pl.LazyFrame({"row_nr": list(range(row_count))})

    expr = pl.lit(False)

    if 'f' in valid_site_types:
        expr |= pl.col("Nonsingleton_Alleles") == 1
    if 'b' in valid_site_types:
        expr |= pl.col("Nonsingleton_Alleles") == 2
    if 't' in valid_site_types:
        expr |= pl.col("Nonsingleton_Alleles") == 3
    if 'q' in valid_site_types:
        expr |= pl.col("Nonsingleton_Alleles") == 4
    if 'p' in valid_site_types:
        expr |= pl.col("Nonsingleton_Alleles") == 5
    if 'u' in valid_site_types:
        expr |= (pl.col("Singletons") > 0) & (pl.col("Nonsingleton_Alleles") == 1)
    if 's' in valid_site_types:
        expr |= (pl.col("Singletons") > 0) & (pl.col("Nonsingleton_Alleles") == 0)

    pl.concat([sites,row_numbers],how='horizontal').filter(expr).select("row_nr").collect(streaming=True).write_parquet(valid_sites_file)
    return valid_sites_file

def find_elbow_index(missing, counts):
    p1 = np.array([missing[0], counts[0]])
    p2 = np.array([missing[-1], counts[-1]])

    line_vec = p2 - p1
    denom = np.sqrt(np.sum(line_vec**2))
    if denom == 0:
        return 0
    line_vec_norm = line_vec / denom

    distances = []
    for i in range(len(missing)):
        p = np.array([missing[i], counts[i]])
        vec = p - p1
        proj_len = np.dot(vec, line_vec_norm)
        proj_point = p1 + proj_len * line_vec_norm
        dist = np.sqrt(np.sum((p - proj_point) ** 2))
        distances.append(dist)

    return int(np.argmax(distances))

def get_max_missing(site_file, pass_gross_exclusion_file, pass_site_type_file, temp_directory):

    sites = pl.scan_parquet(site_file)
    row_count = pq.ParquetFile(site_file).metadata.num_rows
    row_numbers = pl.LazyFrame({"row_nr": list(range(row_count))})
    row_sites = pl.concat([sites,row_numbers],how="horizontal")

    gross_pass = pl.scan_parquet(pass_gross_exclusion_file).select("row_nr")
    site_pass = pl.scan_parquet(pass_site_type_file).select("row_nr")
    intersection = gross_pass.join(site_pass, on="row_nr", how="inner")
    
    filtered = (row_sites.join(intersection, on="row_nr", how="semi"))

    grouped_counts = (
        filtered
        .group_by(["Site_Code", "Missing"])
        .agg(pl.len().alias("count"))
        .sort(["Site_Code", "Missing"])
    ).collect(streaming=True)

    if grouped_counts.height == 0:
        return 0

    elbows = {}
    for site_code, group in grouped_counts.group_by("Site_Code"):
        if group.height < 30:
            continue

        missing_vals = group["Missing"].to_numpy()
        counts = group["count"].to_numpy()

        elbow_index = find_elbow_index(missing_vals, counts)
        elbows[site_code] = missing_vals[elbow_index]

    if elbows:
        return int(np.median(list(elbows.values())))
    else:
        return 0

def filter_missing(site_file,max_missing,temp_directory):
    
    valid_missing_file = os.path.join(temp_directory, "Pass_Missing.parquet")
    sites = pl.scan_parquet(site_file)
    row_count = pq.ParquetFile(site_file).metadata.num_rows
    row_numbers = pl.LazyFrame({"row_nr": list(range(row_count))})
    row_sites = pl.concat([sites,row_numbers],how="horizontal")

    (
        row_sites
        .filter(pl.col('Missing') <= max_missing)
        .select("row_nr")
        .collect(streaming=True)
        .write_parquet(valid_missing_file)
    )

    return valid_missing_file

def merge_filters(pass_gross_exclusion_file,pass_site_type_file,pass_missing_file,temp_directory):

    final_filter_file = os.path.join(temp_directory,"Final_Filter.parquet")
    gross = pl.scan_parquet(pass_gross_exclusion_file).select("row_nr")
    site = pl.scan_parquet(pass_site_type_file).select("row_nr")
    missing = pl.scan_parquet(pass_missing_file).select("row_nr")

    intersect_rows = (
        gross
        .join(site, on="row_nr", how="inner")
        .join(missing, on="row_nr", how="inner")
        .collect(streaming=True)
        .write_parquet(final_filter_file)
    )

    return final_filter_file
    
def save_filtered_parquet(raw_parquet,new_parquet,row_file):
    
    included_rows = pl.scan_parquet(row_file).select("row_nr")

    raw_file = pl.scan_parquet(raw_parquet)
    row_count = pq.ParquetFile(raw_parquet).metadata.num_rows
    row_numbers = pl.LazyFrame({"row_nr": list(range(row_count))})
    row_sites = pl.concat([raw_file,row_numbers],how="horizontal")

    (
        row_sites
        .join(included_rows, on="row_nr", how="inner")
        .drop('row_nr')
        .collect(streaming=True)
        .write_parquet(new_parquet)
    )


def filter_sites(site_file,new_site_file,new_missing_file,final_filter_file):

    included_rows = pl.scan_parquet(final_filter_file).select("row_nr")
    lazy_sites = pl.scan_parquet(site_file)
    row_count = pq.ParquetFile(site_file).metadata.num_rows
    row_numbers = pl.LazyFrame({"row_nr": list(range(row_count))})
    row_sites = pl.concat([lazy_sites,row_numbers],how="horizontal")
    
    (
        row_sites
        .join(included_rows,on='row_nr',how="inner")
        .drop('row_nr')
        .collect(streaming=True)
        .write_parquet(new_site_file)
    )

    code_labels = {
        0: "Other",
        1: "Pure_Fixed",
        2: "Pure_Biallelic",
        3: "Pure_Triallelic",
        4: "Pure_Quadallelic",
        5: "Pure_Pentallelic",
        6: "Fixed_wSingleton",
        7: "Biallelic_wSingleton",
        8: "Triallelic_wSingleton",
        9: "Quadallelic_wSingleton",
        10: "Pentallelic_wSingleton"
    }

    df = ( pl.scan_parquet(new_site_file) .select(['Missing', 'Site_Code']) .group_by(["Missing", "Site_Code"]) .len() .collect() )   
    
    missing_summary = {}
 
    for row in df.iter_rows(named=True):
        m = row["Missing"]
        c = row["Site_Code"]
        v = row["len"]
        col_name = code_labels.get(c, f"Code_{c}")
        missing_summary.setdefault(m, {})[col_name] = v

    wide_df = (
        pl.DataFrame([
            {"Missing": m, **counts}
            for m, counts in missing_summary.items()
        ])
        .fill_null(0)
    )

    nonzero_cols = ["Missing"] + [
        col for col in wide_df.columns if col != "Missing" and wide_df[col].sum() > 0
    ]

    wide_df.select(nonzero_cols).sort("Missing").to_pandas().to_csv(new_missing_file, sep="\t", index=False)

##### MAIN #####

# region 00: Parse args and set up directories
args = parse_args()

filter_id = args.filter_id
join_id = args.join_id

joined_directory = os.path.abspath(args.joined_dir)
output_directory = os.path.abspath(args.output_directory)

ref_fasta = os.path.abspath(args.ref_fasta)

scaffold_file = os.path.join(joined_directory, f"{join_id}_Scaffold.parquet")
code_file = os.path.join(joined_directory, f"{join_id}_Codes.parquet")
site_file = os.path.join(joined_directory, f"{join_id}_Sites.parquet")
base_file = os.path.join(joined_directory, f"{join_id}_Bases.parquet")
missing_file = os.path.join(joined_directory, f"{join_id}_Missing.tsv")

for f in [joined_directory,output_directory,ref_fasta,scaffold_file,code_file,site_file,base_file,missing_file]:
    if not os.path.exists(f):
        sys.exit(f"Expected file/directory {f} does not exist...")

# Filter prefix
new_scaffold_file = os.path.join(output_directory, f"{filter_id}_Scaffold.parquet")
new_code_file = os.path.join(output_directory, f"{filter_id}_Codes.parquet")
new_site_file = os.path.join(output_directory, f"{filter_id}_Sites.parquet")
new_base_file = os.path.join(output_directory, f"{filter_id}_Bases.parquet")
new_missing_file = os.path.join(output_directory, f"{filter_id}_Missing.tsv")
filter_json_file = os.path.join(output_directory, f"{filter_id}.json")

for f in [new_scaffold_file,new_code_file,new_site_file,new_base_file,new_missing_file,filter_json_file]:
    if os.path.exists(f):
        sys.exit(f"File {f} already exists...")

# Create temp directory
temp_directory = os.path.join(output_directory, "SNPRS_Temp")
os.mkdir(temp_directory)

# Process reference FASTA
raw_records = [(rec.id, str(rec.seq)) for rec in SeqIO.parse(ref_fasta, "fasta")]
if not raw_records:
    sys.exit("No contigs found.")

raw_records = natsorted(raw_records, key=lambda x: x[0])
total_sites = sum(len(seq) for _, seq in raw_records)

contig_data = [(rec_id.strip().split()[0],len(seq)) for rec_id, seq in raw_records]
index_key = {index: (contig_id,contig_length) for index,(contig_id, contig_length) in enumerate(contig_data)}
contig_key = {contig_id: (index,contig_length) for index,(contig_id, contig_length) in enumerate(contig_data)}

# Get sample count + IDs
sample_ids = pl.scan_parquet(base_file).collect_schema().names()

# endregion

# region 01: Site inclusion/exclusion

requested_site_types = sorted(set(re.sub(r'[^A-Za-z]', '', args.site_types).lower()))

if len(requested_site_types) < 1:
    sys.exit("Error: no valid site types provided.")

invalid = set(requested_site_types) - set("fbtqpsu")

if invalid:
    sys.exit(f"Error: invalid site type(s) [Allowed: fbtqpsu]: {', '.join(sorted(invalid))}")
    
valid_site_types = "".join(requested_site_types)

if (args.no_singletons) and ("s" in valid_site_types or "u" in valid_site_types):
    remove_singletons = False 
else:
    remove_singletons = args.no_singletons

remove_gaps = args.remove_gaps
include_hets = args.include_hets
include_invalid = args.include_invalid

# Missing data
if not args.missing:
    missing_mode = "estimate"
    max_missing = np.nan
elif int(args.missing) == 0:
    max_missing = 0
elif float(args.missing) < 1:
    missing_mode = "prop"
    min_present = int(float(args.missing) * sample_count)
    max_missing = sample_count - min_present
elif int(args.missing) >= 1:
    missing_mode = "value"
    max_missing = int(args.missing)

# endregion

# region 02: Site Filtering

pass_gross_exclusion_file = gross_exclusion(site_file,remove_gaps,include_hets,remove_singletons,include_invalid,temp_directory)
if not os.path.exists(pass_gross_exclusion_file):
    sys.exit("No sites remain after processing gaps, hets, singletons, and invalid data")

pass_gross_exclusion_count = pq.ParquetFile(pass_gross_exclusion_file).metadata.num_rows
if pass_gross_exclusion_count == 0:
    sys.exit("No sites remain after processing gaps, hets, singletons, and invalid data")

pass_site_type_file = filter_site_types(site_file,valid_site_types,temp_directory)
if not os.path.exists(pass_site_type_file):
    sys.exit(f"No sites remain after requesting --site_types {valid_site_types}")

pass_site_type_count = pq.ParquetFile(pass_site_type_file).metadata.num_rows
if pass_site_type_count == 0:
    sys.exit(f"No sites remain after requesting --site_types {valid_site_types}")

# Get sites that pass both filters
if missing_mode == "estimate":
    max_missing = get_max_missing(site_file,pass_gross_exclusion_file,pass_site_type_file,temp_directory)

pass_missing_file = filter_missing(site_file,max_missing,temp_directory)
if not os.path.exists(pass_missing_file):
    sys.exit(f"No sites remain after requesting removing sites with >= {max_missing}")

pass_missing_count = pq.ParquetFile(pass_missing_file).metadata.num_rows
if pass_missing_count == 0:
    sys.exit(f"No sites remain after requesting removing sites with >= {max_missing}")

final_filter_file = merge_filters(pass_gross_exclusion_file,pass_site_type_file,pass_missing_file,temp_directory)
if not os.path.exists(final_filter_file):
    sys.exit(f"No sites in common after applying filters")

pass_final_count = pq.ParquetFile(final_filter_file).metadata.num_rows
if pass_final_count == 0:
    sys.exit(f"No sites in common after applying filters")

# endregion

# region 03: Save filtered output 
save_filtered_parquet(scaffold_file,new_scaffold_file,final_filter_file)
save_filtered_parquet(code_file,new_code_file,final_filter_file)
save_filtered_parquet(base_file,new_base_file,final_filter_file)

# Save sites and regnerate missing output
filter_sites(site_file,new_site_file,new_missing_file,final_filter_file)

# endregion

# region 04: Save Summary
filtering_info = {
    "Joined_Directory":joined_directory,
    "Filtered_Directory":output_directory,
    "Ref_FASTA": ref_fasta,
    "Sample_IDs":sample_ids,
    "Site_Types":valid_site_types,
    "Gaps_Removed":str(remove_gaps),
    "Singletons_Removed":str(remove_singletons),
    "Hets_Included":str(include_hets),
    "Invalid_Included":str(include_invalid),
    "Missing_Arg":str(args.missing),
    "Max_Missing":str(max_missing),
    "Final_Site_Count":str(pass_final_count)
}

with open(filter_json_file, "w", encoding="utf-8") as f:
    json.dump(filtering_info, f, indent=4)

shutil.rmtree(temp_directory)

# endregion