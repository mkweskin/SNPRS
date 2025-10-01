import os
import polars as pl
import argparse
from itertools import combinations, permutations
import pandas as pd
import pyarrow.parquet as pq
import pyarrow as pa
import json
import sys

# Site Types
# 0: Fixed Base
# 1: Fixed Deletion
# 2: Fixed Insertion
# 3: Het with only ACTG
# 4: Het with any -
# 5: Het Insertion
# 6: Ploidy Fail

def parse_args():
    
    parser = argparse.ArgumentParser(description="Based on min_depth, min_freq, and max_alleles, saves a parquet file with called bases.")

    parser.add_argument('-p',dest="sample_parquet", required=True, type=str, default=None, help="Path to sample parquet file")
    parser.add_argument('-o',dest="output_parquet", type=str, default=None, help="Path to output parquet file with called bases")
    parser.add_argument('-min_depth',dest="site_cov", type=int,default=2, help="Minimum total read depth to consider a site [Default: 2]")
    parser.add_argument('-min_support',dest="allele_cov", type=int,default=2, help="Minimum read depth to consider an allele [Default: 2]")
    parser.add_argument('-min_freq',dest="min_freq", type=float,default=0.15, help="Minimum allele frequency to consider an alternative allele [Default: 0.15]")
    parser.add_argument('-max_alleles',dest="max_alleles", type=int,default=2, help="Maximum number of alleles above min_freq allowed [Default: 2 (diploid)]")

    return parser.parse_args()

def make_degen_dict():
    iupac_map = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "AC": "M",
        "AG": "R",
        "AT": "W",
        "CG": "S",
        "CT": "Y",
        "GT": "K",
        "ACG": "V",
        "ACT": "H",
        "AGT": "D",
        "CGT": "B",
        "ACGT": "N",
    }

    bases = ["A", "C", "G", "T"]
    degen_dict = {}

    for n in range(1, len(bases)+1):
        for combo in combinations(bases, n):
            sorted_key = "".join(sorted(combo))

            code = iupac_map.get(sorted_key, None)

            for perm in permutations(combo):
                key = "".join(perm)
                degen_dict[key] = code

    return degen_dict

#### MAIN ####
args = parse_args()

sample_parquet = os.path.abspath(args.sample_parquet)
read_cov = int(args.site_cov)
allele_cov = int(args.allele_cov)
read_cov = max(read_cov, allele_cov)

min_freq = float(args.min_freq)
max_alleles = int(args.max_alleles)

degen_dict = make_degen_dict()

output_parquet = os.path.abspath(args.output_parquet) if args.output_parquet else f"{os.path.splitext(sample_parquet)[0]}_Called.parquet"
if os.path.exists(output_parquet):
    sys.exit(f"{output_parquet} already exists...")
    
# Extract metadata
schema = pq.read_schema(sample_parquet)
metadata_bytes = schema.metadata or {}
og_metadata = {k.decode("utf-8"): v.decode("utf-8") for k, v in metadata_bytes.items()}
new_metadata = og_metadata.copy()
sample_name = og_metadata['sample_id']

# Add QC parameters
new_metadata["min_read_coverage"] = str(read_cov)
new_metadata["min_allele_coverage"] = str(allele_cov)
new_metadata["min_allele_frequency"] = str(min_freq)
new_metadata["ploidy"] = str(max_alleles)

# Add type map
type_code_map = {
    0: "Fixed_Base",
    1: "Fixed_Deletion",
    2: "Fixed_Insertion",
    3: "Het_Base",
    4: "Het_Deletion",
    5: "Het_Insertion"
}

new_metadata["Type_Code_Map"] = "0: Fixed_Base; 1: Fixed_Deletion; 2: Fixed_Insertion; 3: Het_Base; 4: Het_Deletion; 5: Het_Insertion"

base_df_lazy = (
    pl.read_parquet(sample_parquet).lazy()
    .select(["contig_index", "contig_position", "depth", "base", "frequency"])
    .with_columns([
        (pl.col("depth") >= read_cov).alias("pass_depth"),
        (pl.col("frequency") >= min_freq).alias("pass_freq"),
        ((pl.col("depth") * pl.col("frequency")).floor() >= allele_cov).alias("pass_allele"),
    ])
    .with_columns(
        pl.when(~pl.col("pass_depth"))
            .then(pl.lit("Fail_Depth"))
        .when(~pl.col("pass_freq") | ~pl.col("pass_allele"))
            .then(pl.lit("Fail_Allele"))
        .when(pl.col("base").str.starts_with("N"))
            .then(pl.lit("Read_N"))
        .when(pl.col("base").str.starts_with("+"))
            .then(pl.lit("Insertion"))
        .when(pl.col("base").str.starts_with("-"))
            .then(pl.lit("Deletion"))
        .otherwise(pl.lit("Pass"))
        .alias("status")
    )
)

ploidy_fail_positions = (
    base_df_lazy
    .filter(
        (pl.col("status").is_in(["Pass", "Deletion"]))
    )
    .group_by(["contig_index", "contig_position"])
    .agg(pl.len().alias("allele_count"))
    .filter(pl.col("allele_count") > max_alleles)
    .select(["contig_index", "contig_position"])
)

base_df = (
    base_df_lazy.join(
        ploidy_fail_positions.with_columns(pl.lit("Fail_Ploidy").alias("ploidy_status")),
        on=["contig_index", "contig_position"],
        how="left"
    )
    .with_columns(
        pl.when(pl.col("ploidy_status").is_not_null())
        .then(pl.col("ploidy_status"))
        .otherwise(pl.col("status"))
        .alias("status")
    )
    .drop("ploidy_status")
).collect()

status_counts = {
    "Fail_Depth": (
        base_df
        .filter(pl.col("status") == "Fail_Depth")
        .select(["contig_index", "contig_position"])
        .unique()
        .height
    ),
    "Fail_Ploidy": (
        base_df
        .filter(pl.col("status") == "Fail_Ploidy")
        .select(["contig_index", "contig_position"])
        .unique()
        .height
    ),
    **{
        status: (
            base_df
            .filter(pl.col("status") == status)
            .height
        )
        for status in ["Fail_Allele", "Read_N","Insertion","Deletion","Pass"]
    }
}

status_counts_list = [f"{k}:{v}" for k, v in status_counts.items()]
status_count_string = ", ".join(status_counts_list)
new_metadata['Raw_Status_Counts'] = status_count_string

# Process insertions
insertion_df = (
    base_df
    .filter(pl.col("status") == "Insertion")
    .rename({"base": "final_base"})
    .with_columns(
        pl.when(pl.col("frequency") >= 1 - min_freq)
        .then(2)
        .otherwise(5)
        .alias("type")
    )
    .select(["contig_index", "contig_position", "final_base", "type"])
).cast({
        "contig_index": pl.Int64,
        "contig_position": pl.Int64,
        "final_base": pl.Utf8,
        "type": pl.Int64
      })

# Process pass QC alleles
pass_df = (
    base_df
    .filter(
        ((pl.col("status") == "Pass") |
         (pl.col("status") == "Deletion"))
        )
)

pass_counts = (
    pass_df
    .group_by(["contig_index", "contig_position"])
    .agg(pl.len().alias("row_count"))
)

# Get fixed sites
fixed_df = (
    pass_counts
    .filter(pl.col("row_count") == 1)
    .select(["contig_index", "contig_position"])
    .join(pass_df, on=["contig_index", "contig_position"], how="left")
    .rename({"base": "final_base"})
    .select(["contig_index", "contig_position", "final_base"])
    .with_columns(
        pl.when(pl.col("final_base") == "-")
        .then(pl.lit(1))
        .otherwise(pl.lit(0))
        .alias("type")
    )
).cast({
        "contig_index": pl.Int64,
        "contig_position": pl.Int64,
        "final_base": pl.Utf8,
        "type": pl.Int64
      })

# Get het sites
het_df = (
    pass_counts
    .filter(pl.col("row_count") > 1)
    .join(pass_df, on=["contig_index", "contig_position"], how="left")
)


if het_df.height > 0:
    
    base_only_df = (
        het_df
        .filter(pl.col("base").is_in(["A", "C", "G", "T"]))
        .select(['contig_index','contig_position','base'])
        .group_by(["contig_index", "contig_position"])
        .agg(pl.col('base')).with_columns(pl.col('base').list.join('').alias('base_string'))
        .with_columns([ pl.col("base_string").replace_strict(degen_dict).alias("degen_code") ])
        )
    
    gap_positions = (
        het_df
        .filter(pl.col("base").str.contains(r"[-]"))
        .select(["contig_index","contig_position"])
    )
    
    het_base_df = (
        base_only_df
        .filter(~pl.struct(["contig_index", "contig_position"]).is_in(gap_positions))
        .with_columns([
            pl.col("degen_code").alias("final_base"),
            pl.lit(3).alias("type")
        ])
    )
    
    het_gap_df = (
        base_only_df
        .filter(pl.struct(["contig_index", "contig_position"]).is_in(gap_positions))
        .with_columns([
            pl.col("degen_code").str.to_lowercase().alias("final_base"),
            pl.lit(4).alias("type")
        ])
    )
        
    het_df = pl.concat([het_base_df,het_gap_df]).select(['contig_index','contig_position','final_base','type']).cast({
        "contig_index": pl.Int64,
        "contig_position": pl.Int64,
        "final_base": pl.Utf8,
        "type": pl.Int64
      })
    
else:
    het_df = pl.DataFrame(schema={
        "contig_index": pl.Int64,
        "contig_position": pl.Int64,
        "final_base": pl.Utf8,
        "type": pl.Int64
    })

combined_df = pl.concat([insertion_df, fixed_df, het_df]).sort(['contig_index','contig_position'])

# Get final counts
map_df = pl.DataFrame({
    "type": list(type_code_map.keys()),
    "type_label": list(type_code_map.values())
})

type_counts = (
    combined_df.group_by("type").agg(pl.len().alias("count"))
    .join(map_df, on="type", how="left")
    .sort("type")
)

type_counts_dict = {
    (row["type_label"]): row["count"]
    for row in type_counts.select(["type_label", "count"]).to_dicts()
}

type_counts_list = [f"{k}:{v}" for k, v in type_counts_dict.items()]
type_counts_str = ", ".join(type_counts_list)
new_metadata['Called_Base_Counts'] = type_counts_str

# Save final calls
final_metadata = {k.encode(): v.encode() for k, v in new_metadata.items()}
final_arrow = combined_df.to_arrow().replace_schema_metadata(final_metadata)
pq.write_table(final_arrow, output_parquet, compression="snappy")
