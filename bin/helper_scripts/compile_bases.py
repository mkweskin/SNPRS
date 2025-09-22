import polars as pl
from natsort import natsorted
import sys
import os
import re
import pyarrow.parquet as pq
import pyarrow as pa
import json

def compileBases(base_parquet,base_file):

    writer = None
        
    with open(base_file, "r") as f:
        chunk_files = [line.strip() for line in f if line.strip()]
    
    for chunk_file in chunk_files:
        
        with open(chunk_file, "r") as f:
            chunk_list = json.load(f)

        sample_id_paths = [(os.path.basename(path).split("_Called")[0], path) for path in chunk_list]
        sorted_path_sample_pairs = natsorted(sample_id_paths, key=lambda x: x[0])

        dfs = [
            pl.scan_parquet(path)
            .select("final_base")
            .rename({"final_base": sample_id})
            .fill_null("N") # May want to keep as null?
            for sample_id, path in sorted_path_sample_pairs
        ]

        df = pl.concat(dfs, how="horizontal").collect().to_arrow()

        if writer is None:
            writer = pq.ParquetWriter(base_parquet, df.schema, compression="snappy")

        writer.write_table(df)

    if writer is not None:
        writer.close()
        
if __name__ == "__main__":
    base_parquet,base_file = sys.argv[1],sys.argv[2]
    compileBases(base_parquet,base_file)