import argparse
import polars as pl
import sys

parser = argparse.ArgumentParser()

parser.add_argument("--in",dest="scaffold_parquet",type=str,required=True,help="Path to scaffold base parquet")
parser.add_argument("--out",dest="output_parquet",type=str,required=True,help="Path to filtered output parquet")

parser.add_argument("--alleles",dest="allele_list",type=str,default="2,3,4,5",help="Comma-separated list of 1-5 (1=fixed, 2=biallelic, 3=triallelic, etc.) [Default: 2,3,4,5]")

parser.add_argument("--min_fixed",dest="min_fixed",type=int,default=None,help="Required number of samples with fixed nucleotides [Default: Any]")
parser.add_argument("--max_ploidy",dest="max_ploidy",type=int,default=None,help="Maximum number of sites allowed with ploidy failure [Default: Any]")
parser.add_argument("--min_clade",dest="min_clade",type=int,default=None,help="Only return SNPs where the second most common allele occurs at least this many times [Default: Any]")

parser.add_argument("--no_sing",action="store_true", help="Remove sites if any sample has singleton variation [Default: False]")
parser.add_argument("--no_gap",action="store_true",help="Remove sites if any sample has an indel [Default: False]")

args = parser.parse_args()

try:
    allele_list = [int(x.strip()) for x in args.allele_list.split(",")]
except ValueError:
    sys.exit("--alleles must be a comma-separated list of integers")

invalid = set(allele_list) - {1, 2, 3, 4, 5}
if invalid:
    sys.exit("--alleles values must be 1, 2, 3, 4, or 5")

lf = pl.scan_parquet(args.scaffold_parquet)

filter_expr = pl.lit(True)

if args.no_gap:
    filter_expr &= pl.col("gap") == 0
    value_cols = ["a", "c", "t", "g"]
else:
    value_cols = ["a", "c", "t", "g", "gap"]

if args.no_sing:
    filter_expr &= ~pl.any_horizontal(
        [pl.col(col) == 1 for col in value_cols]
    )

if args.min_fixed is not None:
    filter_expr &= pl.col("fixed") >= args.min_fixed

if args.max_ploidy is not None:
    filter_expr &= pl.col("pf") <= args.max_ploidy

if args.min_clade is not None:
    max2 = (
        pl.concat_list(value_cols)
        .list.sort(descending=True)
        .list.get(1)
    )
    filter_expr &= max2 >= args.min_clade

filter_expr &= pl.col("pi_alleles").is_in(allele_list)

(
    lf.filter(filter_expr)
      .sort(["contig_index", "contig_position"])
      .sink_parquet(
          args.output_parquet,
          compression="snappy",
      )
)