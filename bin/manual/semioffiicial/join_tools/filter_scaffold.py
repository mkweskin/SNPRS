import argparse
import polars as pl
import sys
from pathlib import Path
import os
import math

parser = argparse.ArgumentParser()

parser.add_argument("--in",dest="scaffold_parquet",type=str,required=True,help="Path to scaffold base parquet")
parser.add_argument("--filter_id",dest="filter_id",type=str,required=True,help="Name for filtered dataset")

parser.add_argument("--alleles",dest="allele_list",type=str,default="2,3,4,5",help="Comma-separated list of 1-5 (1=fixed, 2=biallelic, 3=triallelic, etc.) [Default: 2,3,4,5]")

parser.add_argument("--covered",dest="min_covered",type=float,default=None,help="Required number of samples with coverage [Default: Any]")
parser.add_argument("--fixed",dest="min_fixed",type=float,default=None,help="Required number of samples with fixed nucleotides [Default: Any]")
parser.add_argument("--het",dest="min_het",type=float,default=None,help="Required number of samples with heterozygous nucleotides [Default: Any]")
parser.add_argument("--ploidy",dest="max_ploidy",type=float,default=None,help="Maximum number of sites allowed with ploidy failure [Default: Any]")

parser.add_argument("--min_clade",dest="min_clade",type=int,default=None,help="Only return SNPs where the second most common allele occurs at least this many times [Default: Any]")

parser.add_argument("--no_sing",action="store_true", help="Remove sites if any sample has singleton variation [Default: False]")
parser.add_argument("--no_gap",action="store_true",help="Remove sites if any sample has an indel [Default: False]")
parser.add_argument("--no_het",action="store_true",help="Remove sites if any sample has a heterozygous call [Default: False]")

args = parser.parse_args()

output_directory = Path(args.scaffold_parquet).parent
filter_id = args.filter_id

filter_file = os.path.join(output_directory,f"{filter_id}_Filtered.parquet")
if os.path.exists(filter_file):
    sys.exit(f"{filter_file} exists...")

lf = pl.scan_parquet(args.scaffold_parquet)

sample_count = 49851
#sample_count = (
#    lf
#    .select((pl.col("cov") + pl.col("uncovered")).alias("total"))
#    .head(1)
#   .collect(engine="streaming")
#   .item()
#)

filter_expr = pl.lit(True)

base_cols = ["a", "c", "t", "g", "gap"]

if args.no_het:
    filter_expr &= pl.col("het") == 0
    
if args.no_gap:
    filter_expr &= pl.col("gap") == 0

if args.no_sing:
    filter_expr &= ~pl.any_horizontal(
        [pl.col(col) == 1 for col in base_cols]
    )

try:
    allele_list = [int(x.strip()) for x in args.allele_list.split(",")]
except ValueError:
    sys.exit("--alleles must be a comma-separated list of integers")
    
invalid = set(allele_list) - {1, 2, 3, 4, 5}
if invalid:
    sys.exit("--alleles values must be 1, 2, 3, 4, or 5")

filter_expr &= pl.col("pi_alleles").is_in(allele_list)

if args.min_covered is not None:
    if 0 < args.min_covered < 1:
        cov_threshold = math.floor(sample_count * args.min_covered)
    else:
        cov_threshold = int(args.min_covered)
    filter_expr &= pl.col("cov") >= cov_threshold

if args.min_fixed is not None:
    if 0 < args.min_fixed < 1:
        fixed_threshold = math.floor(sample_count * args.min_fixed)
    else:
        fixed_threshold = int(args.min_fixed)

    filter_expr &= pl.col("fixed") >= fixed_threshold

if args.min_het is not None:
    if 0 < args.min_het < 1:
        het_threshold = math.floor(sample_count * args.min_het)
    else:
        het_threshold = int(args.min_het)

    filter_expr &= pl.col("het") >= het_threshold

if args.max_ploidy is not None:
    if 0 < args.max_ploidy < 1:
        ploidy_threshold = math.ceil(sample_count * args.max_ploidy)
    else:
        ploidy_threshold = int(args.max_ploidy)

    filter_expr &= pl.col("pf") < ploidy_threshold

if args.min_clade is not None:
    max2 = (
        pl.concat_list(base_cols)
        .list.sort(descending=True)
        .list.get(1)
    )
    filter_expr &= max2 >= args.min_clade

(
    lf.filter(filter_expr)
      .sort(["contig_index", "contig_position"])
      .sink_parquet(
          filter_file,
          compression="snappy",
      )
)

print(filter_file)