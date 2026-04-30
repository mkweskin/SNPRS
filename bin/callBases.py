import os
import polars as pl
import argparse
from itertools import combinations, permutations
import pandas as pd
import pyarrow.parquet as pq
import pyarrow as pa
import json
import sys

def parse_args():
    
    parser = argparse.ArgumentParser(description="Based on min_depth, min_freq, and max_alleles, saves a parquet file with called bases.")

    parser.add_argument('-p',dest="sample_parquet", required=True, type=str, default=None, help="Path to sample parquet file")
    parser.add_argument('-o',dest="output_parquet", type=str, default=None, help="Path to output parquet file with called bases")
    parser.add_argument('-min_depth',dest="site_cov", type=int,default=1, help="Minimum total read depth to consider a site [Default: 1]")
    parser.add_argument('-min_support',dest="allele_cov", type=int,default=1, help="Minimum read depth to consider an allele [Default: 1]")
    parser.add_argument('-min_freq',dest="min_freq", type=float,default=0.15, help="Minimum allele frequency to consider an alternative allele [Default: 0.15]")
    parser.add_argument('-max_alleles',dest="max_alleles", type=int,default=2, help="Maximum number of alleles above min_freq allowed [Default: 2 (diploid)]")
    return parser.parse_args()

def make_base_int_map():

    bases = ["A", "C", "G", "T"]
    gap = "-"

    base_int = {}
    int_base = {}
    counter = 1

    for r in range(1, len(bases) + 1):
        for combo in combinations(bases, r):
            key = "".join(sorted(combo))
            base_int[key] = counter
            int_base[counter] = set(combo)
            counter += 1

    base_int[gap] = counter
    int_base[counter] = {gap}
    counter += 1

    for r in range(1, len(bases) + 1):
        for combo in combinations(bases, r):
            baseset = set(combo) | {gap}
            key = "".join(sorted(baseset))
            base_int[key] = counter
            int_base[counter] = baseset
            counter += 1

    return base_int, int_base

#### MAIN ####
args = parse_args()

sample_parquet = os.path.abspath(args.sample_parquet)
read_cov = int(args.site_cov)
allele_cov = int(args.allele_cov)

min_freq = float(args.min_freq)
max_alleles = int(args.max_alleles)

output_parquet = os.path.abspath(args.output_parquet) if args.output_parquet else f"{os.path.splitext(sample_parquet)[0]}_Called.parquet"
if os.path.exists(output_parquet):
    sys.exit(f"{output_parquet} already exists...")
    
# Extract metadata
schema = pq.read_schema(sample_parquet)
metadata_bytes = schema.metadata or {}
og_metadata = {k.decode("utf-8"): v.decode("utf-8") for k, v in metadata_bytes.items()}
new_metadata = og_metadata.copy()
sample_name = og_metadata['sample_id']

# Create base dictionaries
base_int,int_base = make_base_int_map()

# Scan bases
raw_base_df = (
    pl.read_parquet(sample_parquet).lazy()
    .with_columns(
        pl.when(pl.col("depth") < read_cov)
            .then(pl.lit("Fail_Depth"))
        .when(pl.col("frequency") < min_freq)
            .then(pl.lit("Fail_Frequency"))
        .when((pl.col("depth") * pl.col("frequency")).floor() < allele_cov)
            .then(pl.lit("Fail_Allele_Cov"))
        .when(pl.col("base").str.starts_with("N"))
            .then(pl.lit("Read_N"))
        .when(
            (pl.col("base").str.starts_with("+")) &
            (pl.col("frequency") > 1 - min_freq)
        )
            .then(pl.lit("Fixed_Insertion"))
        .when(
            (pl.col("base").str.starts_with("+")) &
            (pl.col("frequency") <= 1 - min_freq)
        )
            .then(pl.lit("Het_Insertion"))
        .otherwise(pl.lit("Pass"))
        .alias("status")
    )
    .select(["contig_index", "contig_position", "depth", "base", "frequency", "status"])
).collect(engine="streaming")

# Get sites that failed QC or were insertions
non_base_df = raw_base_df.filter(pl.col("status") != "Pass")

# Get bases and deletions
base_df = raw_base_df.filter(pl.col("status") == "Pass")

# Count valid bases per position
base_counts = (
    base_df
    .group_by(["contig_index", "contig_position"])
    .agg(pl.len().alias("row_count"))
)

# Process fixed sites (single base at contig/pos)
fixed_df = (
    base_counts
    .filter(pl.col("row_count") == 1)
    .join(base_df, on=["contig_index", "contig_position"], how="left")
    .with_columns([
        pl.col("base").replace_strict(base_int).alias("base_code")
    ])
    .select(["contig_index", "contig_position", "base_code"])
).cast({ "contig_index": pl.Int32, "contig_position": pl.Int32, "base_code": pl.Int8 })

# Processed fixed sites (single base at contig/pos)
het_df = (
    base_counts
    .filter(pl.col("row_count") > 1)
    .join(base_df, on=["contig_index", "contig_position"], how="left")
)

if het_df.height == 0:

    het_df = pl.DataFrame(schema={
        "contig_index": pl.Int32,
        "contig_position": pl.Int32,
        "base_code": pl.Int8,
    })

    ploidy_fail_df = het_df.clone()

else:

    het_check = (
        het_df
        .group_by(["contig_index", "contig_position","row_count"])
        .agg(pl.col("base").alias("bases"))
        .with_columns(pl.col("bases").list.sort().list.join('').alias("base_string"))
        .select("contig_index", "contig_position", "row_count", "base_string")
    )

    bad = het_check.filter(~pl.col("base_string").is_in(list(base_int.keys())))

    if bad.height > 0:
        print("Offending base_string values:")
        print(bad)
        sys.exit("replace_strict will fail")

    het_expanded = (
        het_check
        .with_columns(pl.col("base_string").replace_strict(base_int).alias("base_code"))
        .select(["contig_index", "contig_position", "base_code", "row_count"])
    )

    if max_alleles == 5:

        het_df = (
            het_expanded
            .select(["contig_index","contig_position","base_code"])
        ).cast({ "contig_index": pl.Int32, "contig_position": pl.Int32, "base_code": pl.Int8 })


        ploidy_fail_df = pl.DataFrame(schema={
            "contig_index": pl.Int32,
            "contig_position": pl.Int32,
            "base_code": pl.Int8,
        })

    else:

        het_df = (
            het_expanded
            .filter(pl.col("row_count") <= max_alleles)
            .select(["contig_index","contig_position","base_code"])
        ).cast({ "contig_index": pl.Int32, "contig_position": pl.Int32, "base_code": pl.Int8 })
        
        ploidy_fail_df = (
            het_expanded
            .filter(pl.col("row_count") > max_alleles)
            .select(["contig_index", "contig_position", "base_code"]) 
            .cast({ "contig_index": pl.Int32, "contig_position": pl.Int32, "base_code": pl.Int8 })
            .with_columns([(pl.col("base_code") * -1).alias("base_code")])
        )

# Add QC parameters to metadata
new_metadata["min_read_coverage"] = str(read_cov)
new_metadata["min_allele_coverage"] = str(allele_cov)
new_metadata["min_allele_frequency"] = str(min_freq)
new_metadata["ploidy"] = str(max_alleles)

# Add final counts to metadata
new_metadata["Fail_Depth"] = str(int(non_base_df.filter(pl.col("status") == "Fail_Depth").select(["contig_index", "contig_position"]).unique().height))
new_metadata["Fail_Frequency"] = str(int(non_base_df.filter(pl.col("status") == "Fail_Frequency").height))
new_metadata["Fail_Allele_Cov"] = str(int(non_base_df.filter(pl.col("status") == "Fail_Allele_Cov").height))
new_metadata["Read_N"] = str(int(non_base_df.filter(pl.col("status") == "Read_N").height))
new_metadata["Fixed_Insertion"] = str(int(non_base_df.filter(pl.col("status") == "Fixed_Insertion").height))
new_metadata["Het_Insertion"] = str(int(non_base_df.filter(pl.col("status") == "Het_Insertion").height))
new_metadata['Fixed_Sites'] = str(int(fixed_df.height))
new_metadata['Valid_Het_Sites'] = str(int(het_df.height))
new_metadata['Ploidy_Fail_Sites'] = str(int(ploidy_fail_df.height))

# Save final calls
final_metadata = {k.encode(): v.encode() for k, v in new_metadata.items()}
arrow = (
    pl.concat([fixed_df, het_df, ploidy_fail_df])
      .sort(["contig_index", "contig_position"])
      .to_arrow()
)
arrow = arrow.cast(arrow.schema.with_metadata(final_metadata))
pq.write_table(arrow, output_parquet, compression="snappy")