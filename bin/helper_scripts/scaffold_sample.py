import sys
import polars as pl
import os
from natsort import natsorted
import pyarrow.parquet as pq
import argparse
import shutil
import subprocess
import csv
from multiprocessing import Lock

lock = Lock()

def save_sample_parquet(called_base_parquet, scaffold_parquet, sample_parquet, sample_id,summary_file):

    col = sample_id

    lazy_scaffold = pl.scan_parquet(scaffold_parquet)
    lazy_called = pl.scan_parquet(called_base_parquet)
    lazy_called_sites = lazy_called.select(['contig_index','contig_position'])

    fixed_codes = pl.Series([1, 2, 3, 4, 16])

    missing_df = (
        lazy_scaffold
        .join(lazy_called_sites, on=['contig_index','contig_position'], how='anti')
        .with_columns([pl.lit(0).alias(sample_id)])
        .select(['contig_index','contig_position',sample_id])
        .cast({ "contig_index": pl.Int32, "contig_position": pl.Int32, sample_id: pl.Int8 })
    )

    called_df = (
        lazy_scaffold
        .join(lazy_called, on=['contig_index','contig_position'])
        .select(['contig_index','contig_position','base_code'])
        .rename({"base_code": sample_id})
        .cast({ "contig_index": pl.Int32, "contig_position": pl.Int32, sample_id: pl.Int8 })
    )

    missing_count = missing_df.select(pl.len()).collect().item()
    
    called_count = called_df.filter(pl.col(sample_id) > 0).select(pl.len()).collect().item()
    
    fixed_count = (
        called_df
        .filter(pl.col(sample_id).is_in(fixed_codes))
        .select(pl.len())
        .collect()
        .item()
    )

    ploidy_fail_count = (
        called_df
        .filter(pl.col(sample_id) < 0)
        .select(pl.len())
        .collect()
        .item()
    )

    het_count = called_count - fixed_count

    if os.path.exists(summary_file):
        with lock, open(summary_file, "a", newline="") as out_f:
            writer = csv.writer(out_f, delimiter="\t")
            writer.writerow([
                sample_id,
                fixed_count,
                het_count,
                ploidy_fail_count,
                missing_count
            ])

    (
        pl.concat([called_df, missing_df])
        .sort(['contig_index','contig_position'])
        .select([sample_id])
        .collect(streaming=True)
        .write_parquet(sample_parquet, compression="snappy")
    )

def parse_args():
    parser = argparse.ArgumentParser(description="Based off a scaffold parquet, get base information for a sample")
    
    parser.add_argument("--called_bases", dest="called_base_file", type=str, required=True,help="File with paths to 2+ Called_Bases parquet")
    parser.add_argument("--scaffold", dest="scaffold_file", type=str, required=True,help="Path to scaffold parquet")
    parser.add_argument("--join_id", dest="join_id", type=str, required=True,help="Prefix for output files")
    parser.add_argument("--out_dir", dest="output_directory", type=str, default=None,help="Path to output directory [Default: cwd]")

    return parser.parse_args()

# region 00: Parse args
args = parse_args()

join_id = str(args.join_id)

# Called base file
called_base_file = os.path.abspath(args.called_base_file)
if not os.path.exists(called_base_file):
    sys.exit(f"{called_base_file} (--called_bases) does not exist")

schema = pq.read_schema(called_base_file)
metadata_bytes = schema.metadata or {}
og_metadata = {k.decode("utf-8"): v.decode("utf-8") for k, v in metadata_bytes.items()}

sample_id = og_metadata['sample_id']

# Scaffold file
scaffold_file = os.path.abspath(args.scaffold_file)
if not os.path.exists(scaffold_file):
    sys.exit(f"{scaffold_file} (--scaffold_file) does not exist")

# Output directory
if args.output_directory is None:
    output_directory = os.getcwd()
else:
    output_directory = os.path.abspath(args.output_directory)
if not os.path.exists(output_directory):
    sys.exit(f"{output_directory} (--out_dir) does not exist")

temp_directory = os.path.join(output_directory,f"Temp_{join_id}")
if not os.path.exists(temp_directory):
    os.mkdir(temp_directory)

# Output files
summary_file = os.path.join(output_directory, f"{join_id}_Site_Counts.tsv")
sample_parquet = os.path.join(temp_directory,f"Scaffolded_{sample_id}.parquet")

# endregion

# region 01: Scaffold sample

save_sample_parquet(called_base_file,scaffold_file,sample_parquet,sample_id,summary_file)

# endregion
