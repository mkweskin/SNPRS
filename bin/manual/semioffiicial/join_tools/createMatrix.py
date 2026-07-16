#!/usr/bin/env python3

import argparse
import os
import polars as pl


def parse_args():
    p = argparse.ArgumentParser(description="Create an empty sample-by-site matrix.")
    p.add_argument("-s", "--scaffold", required=True, help="Scaffold parquet file")
    p.add_argument("-c", "--count", type=int, required=True, help="Sample count")
    p.add_argument("-n", "--name", default="all_samples", help="Output prefix")
    return p.parse_args()


def main():
    args = parse_args()

    scaffold_file = os.path.abspath(args.scaffold)
    scaffold_dir = os.path.dirname(scaffold_file)

    output_directory = os.path.join(scaffold_dir, args.name)
    os.makedirs(output_directory, exist_ok=True)

    matrix_file = os.path.join(output_directory, f"{args.name}_Matrix.bin")
    called_base_file = os.path.join(output_directory, f"{args.name}_Called_Bases.txt")
    
    if os.path.exists(matrix_file):
        raise RuntimeError(
            f"{matrix_file} already exists..."
        )

    n_sites = (
        pl.scan_parquet(scaffold_file)
        .select(pl.len())
        .collect()
        .item()
    )

    with open(matrix_file, "wb") as f:
        f.truncate(args.count * n_sites)
    
    print(",".join([
        output_directory.strip(),
        scaffold_file.strip(),
        matrix_file.strip(),
        called_base_file.strip()
    ]))

if __name__ == "__main__":
    main()