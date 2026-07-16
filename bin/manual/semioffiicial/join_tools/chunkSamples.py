#!/usr/bin/env python3

import argparse
import os
import polars as pl
import sys

def parse_args():
    p = argparse.ArgumentParser(description="Combine per-sample BIN files into a single sample-by-site matrix.")
    p.add_argument("-c", "--called_base_file", required=True, help="Path to file with called base paths")
    p.add_argument("-s", "--chunk_size", type=int,default=1, help="Chunk size")
    p.add_argument("-n", "--name", required=True, help="Output prefix")
    return p.parse_args()


def fetch_base_parquets(file_path):

    with open(file_path, "r") as f:
        paths = [os.path.abspath(line.strip())
                 for line in f
                 if line.strip()]

    if not paths:
        sys.exit("Error: You must provide at least one called base parquet file.")

    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        sys.exit("Error: Missing parquet files:\n" + "\n".join("  " + m for m in missing))

    return paths

def main():

    args = parse_args()
    
    called_base_files = fetch_base_parquets(args.called_base_file)
    sample_count = len(called_base_files)
    names = sorted([os.path.basename(f).replace("_Called.parquet","") for f in called_base_files])
    if args.chunk_size > 0:
        chunk_ids = [i // args.chunk_size for i in range(sample_count)]
    else:
        chunk_ids = [0] * sample_count
           
    output_directory = os.path.abspath(os.path.dirname(args.called_base_file))
    names_file = os.path.join(output_directory, f"{args.name}_Chunks.tsv")

    with open(names_file, "w") as out:
        out.write("row_num\tchunk_id\tsample_id\n")
        for i, (chunk, name) in enumerate(zip(chunk_ids, names)):
            out.write(f"{i}\t{chunk}\t{name}\n")
    
    print(names_file,end="")

if __name__ == "__main__":
    main()