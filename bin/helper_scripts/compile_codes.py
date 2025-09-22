import polars as pl
from natsort import natsorted
import sys
import pyarrow.parquet as pq


def compileCodes(output_path, input_list_file):
    with open(input_list_file, "r") as f:
        paths = [line.strip() for line in f if line.strip()]

    sample_paths = natsorted(paths)

    writer = None

    for path in sample_paths:
        df = pl.read_parquet(path).to_arrow()

        if writer is None:
            writer = pq.ParquetWriter(output_path, df.schema, compression="snappy")

        writer.write_table(df)

    if writer is not None:
        writer.close()
        
if __name__ == "__main__":
    output_path, input_list_file = sys.argv[1], sys.argv[2]
    compileCodes(output_path, input_list_file)
    
