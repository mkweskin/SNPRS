#!/usr/bin/env python3
import sys
import polars as pl
import numpy as np
import json
import struct
import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Score chunk of scaffold")    
    parser.add_argument("--called", dest="called_parquet", type=str, required=True, help="Path to called base parquet")
    parser.add_argument("--scaffold", dest="scaffold_parquet", type=str, required=True, help="Path to scaffold base parquet")
    parser.add_argument("--out", dest="out_dir", type=str, help="Path to scaffold base parquet")
    parser.add_argument("--bed", dest="bed_file", type=str, help="4 column tsv with locus information")
    return parser.parse_args()


def code_to_base(x):
    base_convert_dict = {
        1:  "A",
        2:  "C",
        3:  "G",
        4:  "T",
        5:  "M",
        6:  "R",
        7:  "W",
        8:  "S",
        9:  "Y",
        10: "K",
        11: "V",
        12: "H",
        13: "D",
        14: "B",
    }
    return base_convert_dict.get(x, "N")

def wrap80(seq):
    return "\n".join(seq[i:i+80] for i in range(0, len(seq), 80))


args = parse_args()

site_info = pl.scan_parquet(args.scaffold_parquet).select(['contig_index','contig_position'])
called_info = pl.scan_parquet(args.called_parquet).select(['contig_index','contig_position','base_code'])

sample_name = os.path.basename(args.called_parquet)
if sample_name.endswith("_Called.parquet"):
    sample_name = sample_name.replace("_Called.parquet", "")
else:
    sample_name = sample_name.replace(".parquet", "")
    
if args.out_dir:
    output_directory = os.path.abspath(args.out_dir)
else:
    output_directory = f"{os.path.dirname(os.path.abspath(args.scaffold_parquet))}/FASTA"
    
if not os.path.exists(output_directory):
    os.mkdir(output_directory)

joined = (
    site_info
    .join(called_info, on=["contig_index", "contig_position"], how="left")
    .with_columns(
        pl.col("base_code").fill_null(0)
    )
    .sort(['contig_index','contig_position'])
).collect(engine="streaming")

if args.bed_file:
    
    bed_df = pl.read_csv(
        args.bed_file,
        separator="\t",
        has_header=True,
        new_columns=["contig_index", "start", "stop", "name"]
    )

    for locus in bed_df.iter_rows(named=True):

        locus_index = locus["contig_index"]
        locus_start = locus["start"]
        locus_stop = locus["stop"]
        locus_id = locus["name"]
        
        if locus_id:
            locus_df = joined.filter(
                (pl.col("contig_index") == locus_index)
                & (pl.col("contig_position") >= locus_start)
                & (pl.col("contig_position") <= locus_stop)
            ).sort("contig_position")

            bases = [code_to_base(x) for x in locus_df["base_code"]]
            sequence = "".join(bases)

            out_path = f"{output_directory}/{locus_id}_{sample_name}.fasta"

            with open(out_path, "w", newline="\n") as f:
                f.write(f">{sample_name}\n")
                f.write(wrap80(sequence))
                f.write("\n")

else:
    
    bases = [code_to_base(x) for x in joined["base_code"]]
    sequence = "".join(bases)

    out_path = f"{output_directory}/{sample_name}.FASTA"

    with open(out_path, "w", newline="\n") as f:
        f.write(f">{sample_name}\n")
        f.write(wrap80(sequence))
        f.write("\n")
    