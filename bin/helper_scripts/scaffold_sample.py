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

def save_sample_parquet(raw_sample_parquet, called_base_parquet, scaffold_parquet, sample_parquet, sample_id,summary_file):

    valid_sites = [0, 1, 3, 4]
    col = sample_id

    lazy_scaffold = pl.scan_parquet(scaffold_parquet)
    lazy_called = pl.scan_parquet(called_base_parquet).filter(pl.col("type").is_in(valid_sites))

    # Sites in scaffold but not called
    lazy_called_sites = lazy_called.select(['contig_index','contig_position']).unique()
    missing_rows = lazy_scaffold.join(lazy_called_sites, on=['contig_index','contig_position'], how='anti')

    if os.path.exists(raw_sample_parquet):

        # Valid coverage (exclude N and insertions)
        lazy_sample = (
            pl.scan_parquet(raw_sample_parquet)
            .select(['contig_index','contig_position','base'])
            .filter((pl.col("base") != "N") & ~pl.col("base").str.starts_with("+"))
            .select(['contig_index','contig_position'])
            .unique()
        )

        # Uncovered (?): absent in raw
        uncovered_rows = (
            lazy_scaffold
            .join(lazy_sample, on=['contig_index','contig_position'], how='anti')
            .with_columns([pl.lit("?").alias(sample_id),pl.lit(5).alias("type")])
        )

        # Filtered (N): absent in called but present in raw
        filtered_rows = (
            missing_rows
            .join(uncovered_rows, on=['contig_index','contig_position'], how='anti')
            .with_columns([pl.lit("N").alias(sample_id),pl.lit(6).alias("type")])
        )

        missing_df = pl.concat([uncovered_rows, filtered_rows]).with_columns(pl.col("type").cast(pl.Int32))

    else:
        # No raw parquet so all missing are uncovered
        missing_df = missing_rows.with_columns([pl.lit("?").alias(sample_id),pl.lit(5).alias("type")]).with_columns(pl.col("type").cast(pl.Int32))


    # Called bases
    called_rows = (
        lazy_scaffold
        .join(lazy_called, on=['contig_index','contig_position'])
        .select(['contig_index','contig_position','final_base','type'])
        .rename({"final_base": sample_id})
        .with_columns(pl.col("type").cast(pl.Int32))
    )

    result = (
        pl.concat([called_rows, missing_df])
        .sort(['contig_index','contig_position'])
        .collect(streaming=True)
    )

    result.select([sample_id]).write_parquet(sample_parquet, compression="snappy")

    type_counts = (
        result
        .group_by("type")
        .len()
        .to_dict(as_series=False)
    )

    type_dict = dict(zip(type_counts.get("type", []), type_counts.get("len", [])))
    counts = {
        "fixed_base": type_dict.get(0, 0),
        "fixed_gap":  type_dict.get(1, 0),
        "het_base":   type_dict.get(3, 0),
        "het_gap":    type_dict.get(4, 0),
        "uncovered":  type_dict.get(5, 0),
        "filtered":   type_dict.get(6, 0),
    }
    
    if os.path.exists(summary_file):
        with lock, open(summary_file, "a", newline="") as out_f:
            writer = csv.writer(out_f, delimiter="\t")
            writer.writerow([
                sample_id,
                counts["fixed_base"],
                counts["fixed_gap"],
                counts["het_base"],
                counts["het_gap"],
                counts["uncovered"],
                counts["filtered"],
            ])


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
raw_parquet_file = og_metadata['sample_parquet']

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

# Output files
summary_file = os.path.join(output_directory, f"{join_id}_Site_Counts.tsv")
sample_parquet = os.path.join(output_directory,f"Scaffolded_{sample_id}.parquet")

# endregion

# region 01: Assess missing data

save_sample_parquet(raw_parquet_file,called_base_file,scaffold_file,sample_parquet,sample_id,summary_file)

# endregion
