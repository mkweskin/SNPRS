#!/usr/bin/env python3
import sys
import os
import numpy as np
import polars as pl
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Fill in matrix column based on called bases")    
    parser.add_argument("-b", dest="called_parquet", type=str, required=True,help="Path to called base parquet")
    parser.add_argument("-c", dest="chunk_file", type=str, required=True,help="Path to chunk information")
    parser.add_argument("-m", dest="matrix_file", type=str, required=True,help="Path to matrix file")
    parser.add_argument("-i", dest="sample_id", type=str, required=True,help="Name for output file")    
    parser.add_argument("-s", dest="scaffold_file", type=str, required=True,help="Path to scaffold file")    
    return parser.parse_args()

def main():

    args = parse_args()

    matrix_file = os.path.abspath(args.matrix_file)
    called_file = os.path.abspath(args.called_parquet)
    scaffold_file = os.path.abspath(args.scaffold_file)
    
    chunk_df = pl.read_csv(args.chunk_file, separator="\t")
    sample_count = chunk_df.height
        
    sample = (
        chunk_df
        .filter(pl.col("sample_id") == args.sample_id)
        .select(["row_num", "chunk_id"])
    )

    if sample.height != 1:
        raise ValueError(f"Sample '{args.sample_id}' not found in {args.chunk_file}")

    row_number, chunk_id = sample.row(0)
    
    site_info = pl.scan_parquet(scaffold_file).select(['contig_index','contig_position'])
    alignment_length = (
        site_info
        .select(pl.len())
        .collect()
        .item()
    )
    
    called_info = pl.scan_parquet(called_file).select(['contig_index','contig_position','base_code'])

    bases = (
        site_info
        .join(called_info, on=["contig_index", "contig_position"], how="left")
        .with_columns(
            pl.when(pl.col("base_code").is_null())
            .then(0)
            .when(pl.col("base_code").is_in([1, 2, 3, 4, 16]))
            .then(pl.col("base_code"))
            .otherwise(-1)
            .cast(pl.Int8)
            .alias("base_code")
        )
        .sort(['contig_index','contig_position'])
        .select("base_code")
        .collect(streaming=True)
        .get_column("base_code")
        .to_numpy()
    )

    assert len(bases) == alignment_length

    matrix = np.memmap(
        matrix_file,
        dtype=np.int8,
        mode="r+",
        shape=(sample_count, alignment_length)
    )

    matrix[row_number] = bases
    matrix.flush()
    
    print(chunk_id,end="")

if __name__ == "__main__":
    main()
    