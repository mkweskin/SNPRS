import polars as pl
import json
import sys
import os

def save_filtered_join(raw_parquet, new_parquet, row_file):
    
    with open(row_file, "r", encoding="utf-8") as f:
        included_rows = set(json.load(f))

    raw_scan = pl.scan_parquet(raw_parquet).with_row_index("row_nr")
    filtered = raw_scan.filter(pl.col("row_nr").is_in(list(included_rows))).drop("row_nr")
    filtered.collect().write_parquet(new_parquet)
    
if __name__ == "__main__":

    raw_parquet,new_parquet,row_file = sys.argv[1],sys.argv[2],sys.argv[3]
    save_filtered_join(raw_parquet,new_parquet,row_file)
