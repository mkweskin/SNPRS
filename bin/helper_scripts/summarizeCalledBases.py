import pandas as pd
import pyarrow.parquet as pq
import sys
import os
import polars as pl
import argparse
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="Get summary statistics from 1 or more Called.parquet files")

    parser.add_argument("--called_list", dest="called_list", type=str, default = None, help ="Path to a list of _Called.parquet files")
    parser.add_argument("--called_dir", dest="called_dir", type=str,default = None,help="Path to a folder containing _Called.parquet files")
    parser.add_argument("--out",dest="output_file",type=str,required=True,help="Path to output CSV")
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

def parquet_preview(parquet_path,preview_data):

    parquet_file = pq.ParquetFile(parquet_path)
    metadata = parquet_file.schema_arrow.metadata
    
    if metadata is None or len(metadata) == 0:
        print("No metadata found.")
    else:
        print("Metadata:")
        for key, value in metadata.items():
            try:
                print(f"  {key.decode()}: {value.decode()}")
            except Exception:
                print(f"  {key}: (binary data)")

##### MAIN #####

# region 00: Parse args

args = parse_args()

if args.called_list and args.called_dir:
    sys.exit("Cannot run summarizeCalledBases.py with both --called_list and --called_dir")
elif not args.called_list and not args.called_dir:
    sys.exit("Cannot run summarizeCalledBases.py without --called_list or --called_dir")

output_file = os.path.abspath(args.output_file)
if os.path.exists(output_file):
    sys.exit(f"{output_file} already exits...")

if args.called_list:
    called_base_file = os.path.abspath(args.called_list)
    if not os.path.exists(called_base_file):
        sys.exit(f"{called_base_file} (--called_list) does not exist")

    called_base_files = fetch_base_parquets(called_base_file)

else:
    called_dir = os.path.abspath(args.called_dir)
    if not os.path.isdir(called_dir):
        sys.exit(f"{called_dir} (--called_dir) is not a directory")

    called_base_files = sorted(
        glob.glob(os.path.join(called_dir, "*_Called.parquet"))
    )

    if len(called_base_files) == 0:
        sys.exit(f"No *_Called.parquet files found in directory: {called_dir}")


# region 01: Get metadata

metadata_list = []

for called_base in called_base_files:
    parquet_file = pq.ParquetFile(called_base)
    raw_metadata = parquet_file.schema_arrow.metadata or {}

    metadata_dict = {}

    for key, value in raw_metadata.items():
        metadata_dict[key.decode()] = value.decode()

    metadata_list.append(metadata_dict)

df = pd.DataFrame(metadata_list)

def parse_metadata_string(s):

    if not s or not isinstance(s, str):
        return {}
    
    parts = [x.strip() for x in s.split(",")]
    out = {}
    for p in parts:
        if ":" in p:
            key, value = p.split(":", 1)
            out[key.strip()] = value.strip()
    return out

depth_expanded = df["depth_statistics"].apply(parse_metadata_string).apply(pd.Series)
allele_expanded = df["allele_frequencies"].apply(parse_metadata_string).apply(pd.Series)
raw_status_expanded = df["Raw_Status_Counts"].apply(parse_metadata_string).apply(pd.Series)
called_base_expanded = df["Called_Base_Counts"].apply(parse_metadata_string).apply(pd.Series)

final_df = pd.concat(
    [
        df.drop(
            ["depth_statistics", "allele_frequencies", 
             "Raw_Status_Counts", "Called_Base_Counts","Type_Code_Map"],
            axis=1
        ),
        depth_expanded,
        allele_expanded,
        raw_status_expanded,
        called_base_expanded,
    ],
    axis=1
)
final_df.columns = [col.strip().replace(" ", "_").replace("-", "_") for col in final_df.columns]

final_df.to_csv(output_file, sep=",", index=False)