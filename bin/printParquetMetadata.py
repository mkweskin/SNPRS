import pandas as pd
import pyarrow.parquet as pq
import sys
import os
import polars as pl
import argparse

fixed_codes = pl.Series([1, 2, 3, 4, 16, 33,34,35,36,48])

def parse_args():
    parser = argparse.ArgumentParser(description="Print the header metadata from a parquet file")
    parser.add_argument("-p","--p", dest="parquet", type=str, required=True,help="Path to parquet file")
    parser.add_argument("-d","--d", dest="data", action="store_true",help="Preview the top 50 rows of the data table")
    parser.add_argument("--idx", dest="idx",type=int,default = None, help="Preview the top 50 rows of the data table from contig --idx")
    parser.add_argument("--pos", dest="pos",type=int, default = None,help="Preview the top 100 rows of the data table")
    return parser.parse_args()

def parquet_preview(parquet_path,preview_data):

    parquet_file = pq.ParquetFile(parquet_path)
    row_count = parquet_file.metadata.num_rows
    print(f"\nTotal rows: {row_count}")
    
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
            
        if args.idx is not None:
            
            if args.pos is not None:
                idx_df = lf.filter((pl.col("contig_index") == args.idx) & (pl.col("contig_position") == args.pos)).collect(engine="streaming")
            
            else:
                idx_df = lf.filter((pl.col("contig_index") == args.idx)).sort(['contig_index','contig_position']).collect(engine="streaming")
                            
            with pl.Config(tbl_cols=-1, tbl_width_chars=-1,tbl_rows=-1):
                print(idx_df)
        else:
            head_df = lf.limit(25).collect()
            tail_df = lf.tail(25).collect()

            with pl.Config(tbl_cols=-1, tbl_width_chars=-1, tbl_rows=50):
                print("\n=== TOP 25 ===")
                print(head_df)

                print("\n=== BOTTOM 25 ===")
                print(tail_df)


args = parse_args()
parquet_file = os.path.abspath(args.parquet)
preview_data = args.data
parquet_preview(parquet_file,preview_data)