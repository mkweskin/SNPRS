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
import json



def parse_args():
    parser = argparse.ArgumentParser(description="Based off a scaffold parquet, get base information for a sample")
    
    parser.add_argument("--join_id", dest="join_id", type=str, required=True,help="Prefix for output files")
    parser.add_argument("--out_dir", dest="output_directory", type=str, required=True,help="Path to output directory [Default: cwd]")
    return parser.parse_args()

def summarize_sample(sample_id, code_file):
    fixed_codes = {1,2,3,4,16,33,34,35,36,48}
    singleton_codes = {1,2,3,4,16}
    het_codes = set(range(5, 16)) | set(range(17, 32))

    sample_series = pl.scan_parquet(code_file).select(sample_id)

    summary = (
        sample_series
        .with_columns([
            (pl.col(sample_id) == 0).sum().alias("Uncovered"),
            (pl.col(sample_id) < 0).sum().alias("Ploidy_Fail"),
            pl.col(sample_id).is_in(fixed_codes).sum().alias("Fixed"),
            pl.col(sample_id).is_in(het_codes).sum().alias("Hets"),
            pl.col(sample_id).is_in(singleton_codes).sum().alias("Singletons"),
        ])
        .with_columns([
            pl.lit(sample_id).alias("Sample_ID")
        ])
        .select(["Sample_ID","Uncovered","Ploidy_Fail","Fixed","Hets","Singletons"])
        .collect()
    )

    return summary[0]

# region 00: Parse args
args = parse_args()

join_id = args.join_id

join_directory = os.path.abspath(args.output_directory)

if not os.path.exists(join_directory):
    sys.exit(f"{join_directory} (--out_dir) does not exist")

scaffold_parquet = os.path.join(join_directory,f"{join_id}_Scaffold.parquet")
if not os.path.exists(scaffold_parquet):
    sys.exit(f"Expected file {scaffold_parquet} does not exist")

temp_directory = os.path.join(join_directory,f"Temp_{join_id}")
if not os.path.exists(temp_directory):
    sys.exit(f"{temp_directory} does not exist")
    
chunk_tsv = os.path.join(temp_directory,"Chunk_Info.tsv")
if not os.path.exists(chunk_tsv):
    sys.exit(f"Expected file {chunk_tsv} does not exist")

output_code_file = os.path.join(join_directory,f"{join_id}_Codes.parquet")
output_site_file = os.path.join(join_directory,f"{join_id}_Sites.parquet")
output_sample_summary = os.path.join(join_directory,f"{join_id}_Sample_Summary.tsv")
output_missing_summary = os.path.join(join_directory,f"{join_id}_Missing_Summary.tsv")
output_json_file = os.path.join(join_directory,f"{join_id}.json")

if os.path.exists(output_code_file):
    sys.exit(f"{output_code_file} already exists...")

if os.path.exists(output_site_file):
    sys.exit(f"{output_site_file} already exists...")
    
if os.path.exists(output_sample_summary):
    sys.exit(f"{output_sample_summary} already exists...")

if os.path.exists(output_missing_summary):
    sys.exit(f"{output_missing_summary} already exists...")
    
if os.path.exists(output_json_file):
    sys.exit(f"{output_json_file} already exists...")
    
# endregion

# region 01: Compile chunks

try:
    chunk_df = pd.read_csv(chunk_tsv, sep="\t")

    code_files = [
        os.path.join(row["Chunk_Directory"], f"{row['Chunk_ID']}_Codes.parquet")
        for _, row in chunk_df.iterrows()
    ]

    site_files = [
        os.path.join(row["Chunk_Directory"], f"{row['Chunk_ID']}_Sites.parquet")
        for _, row in chunk_df.iterrows()
    ]

    for f in code_files + site_files:
        if not os.path.exists(f):
            sys.exit(f"Expected file {f} does not exist")

    lazy_code_frames = [pl.scan_parquet(f) for f in code_files]
    lazy_site_frames = [pl.scan_parquet(f) for f in site_files]

    pl.concat(lazy_code_frames).sink_parquet(output_code_file, compression="snappy")
    pl.concat(lazy_site_frames).sink_parquet(output_site_file, compression="snappy")

    sample_ids = pl.scan_parquet(output_code_file).collect_schema().names()[2:]

    # endregion

    ### PAUSE ###
    
#    # region 02: Summarize sample data
#    lazy_codes = pl.scan_parquet(output_code_file)
#    sample_ids = lazy_codes.collect_schema().names()[2:]
#
#    sample_summary_list = [summarize_sample(s, output_code_file) for s in sample_ids]
#    pl.concat(sample_summary_list, how="vertical").to_pandas().to_csv(output_sample_summary, sep="\t", index=False)
#
#    # endregion
#
#    # region 03: Summarize site data 
#
#    lazy_sites = pl.scan_parquet(output_site_file)
#
#    code_labels = {
#        0: "Other",
#        1: "Pure_Fixed",
#        2: "Pure_Biallelic",
#        3: "Pure_Triallelic",
#        4: "Pure_Quadallelic",
#        5: "Pure_Pentallelic",
#        6: "Fixed_wSingleton",
#        7: "Biallelic_wSingleton",
#        8: "Triallelic_wSingleton",
#        9: "Quadallelic_wSingleton",
#    }
#
#    df_counts = (
#        lazy_sites
#        .filter(pl.col("Site_Code") != 0)
#        .with_columns([
#            pl.col("Site_Code").replace_strict(code_labels).alias("Site_Type")
#        ])
#        .group_by(["Missing", "Site_Type"])
#        .agg(pl.len().alias("count"))
#        .collect()
#        .to_pandas()
#    )
#
#    summary_pivot = df_counts.pivot_table(
#        values="count",
#        index="Missing",
#        columns="Site_Type",
#        aggfunc="first",
#        fill_value=0
#    ).reset_index()
#
#    summary_pivot.to_csv(output_missing_summary, sep="\t", index=False)

    joined_info = {
        "Join_ID":join_id,
        "Joined_Directory":join_directory,
        #"Sample_IDs":",".join(sample_ids),
        "Scaffold_File":scaffold_parquet,
        "Code_File":output_code_file,
        "Site_File":output_site_file
        #,
        #"Sample_Summary_File":output_sample_summary,
        #"Site_Count_File":output_missing_summary
        }

    with open(output_json_file, "w", encoding="utf-8") as f:
        json.dump(joined_info, f, indent=4)
    
    pass

except Exception as e:
    print(f"Exception occurred: {e}")

else:
    shutil.rmtree(temp_directory)
