import pandas as pd
import pyarrow.parquet as pq
import sys
import os
import polars as pl
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Print the header metadata from a parquet file")
    parser.add_argument("-p","--p", dest="parquet", type=str, required=True,help="Path to parquet file")
    parser.add_argument("-d","--d", dest="data", action="store_true",help="Preview the top 100 rows of the data table")
    return parser.parse_args()

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

    if preview_data:
        lf = pl.read_parquet(parquet_path).lazy()

        row_count = lf.select(pl.len()).collect().item() 
        print(f"\nTotal rows: {row_count}")

        head_df = lf.limit(100).collect()
        print("\nFirst 100 rows:")
        pl.Config.set_tbl_rows(100)
        print(head_df)

args = parse_args()
parquet_file = os.path.abspath(args.parquet)
preview_data = args.data
parquet_preview(parquet_file,preview_data)