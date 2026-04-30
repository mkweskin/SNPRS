import polars as pl
import sys
import textwrap
import os

def make_alignment(base_file,alignment_file):

    sample_ids = pl.scan_parquet(base_file).collect_schema().names()

    for sample_id in sample_ids:
        seq = "".join(pl.read_parquet(base_file, columns=[sample_id])[sample_id].to_list())
        
        if os.path.exists(alignment_file):
            with open(alignment_file, "a", encoding="utf-8") as f:
                f.write(f">{sample_id}\n")
                for line in textwrap.wrap(seq, 80):
                    f.write(line + "\n")
        else:
            with open(alignment_file, "w", encoding="utf-8") as f:
                f.write(f">{sample_id}\n")
                for line in textwrap.wrap(seq, 80):
                    f.write(line + "\n")

if __name__ == "__main__":

    base_file,alignment_file = sys.argv[1],sys.argv[2]
    make_alignment(base_file,alignment_file)