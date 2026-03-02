#!/usr/bin/env python3

import polars as pl
import os
import argparse


parser = argparse.ArgumentParser(description='ID called bases from a SNP parquet')
parser.add_argument('-c','--called_bases',dest="called_bases", type=str,required=True, help='Path to a SNPRS called base file or a file with multiple paths.')
parser.add_argument('-s','--snp',dest="snp_parquet", type=str,required=True, help='Path to SNP parquet file')
parser.add_argument('-n','--name',dest="output_file", default=None,type=str, help='Name for output file')

def snpClassifier(snp_df, sample_id, called_base_path):

    fixed_codes = {1, 2, 3, 4, 16}

    degenerate_map = {
        1:  {5,21,6,22,7,23,11,27,12,28,13,29,15,17,31},
        2:  {5,21,8,24,9,25,11,27,12,28,14,30,15,18,31},
        3:  {6,22,8,24,10,26,11,27,13,29,14,30,15,19,31},
        4:  {7,23,9,25,10,26,12,28,13,29,14,30,15,20,31},
        16: {17,18,29,20,21,22,23,24,25,26,27,28,29,30,31},
    }

    def degenerate_match_expr():
        base = pl.col("base_code")
        snp  = pl.col("SNP_Base")
        return (
            pl.when(snp == 1).then(base.is_in(degenerate_map[1]))
            .when(snp == 2).then(base.is_in(degenerate_map[2]))
            .when(snp == 3).then(base.is_in(degenerate_map[3]))
            .when(snp == 4).then(base.is_in(degenerate_map[4]))
            .when(snp == 16).then(base.is_in(degenerate_map[16]))
            .otherwise(False)
        ).alias("Match")

    def safe_prop(numer, denom):
        return (numer / denom).fill_null(0).alias("prop")

    sample_called_bases = pl.read_parquet(called_base_path)
    groups = snp_df["SNP_Group"].unique().to_list()

    sample_rows = []

    for gr in groups:

        gr_snp = snp_df.filter(pl.col("SNP_Group") == gr)

        df = (
            gr_snp.join(sample_called_bases,
                        on=["contig_index", "contig_position"],
                        how="left")
                 .filter(pl.col("base_code").is_not_null())
        )

        fixed_df = df.filter(pl.col("base_code").is_in(fixed_codes)).with_columns(
            (pl.col("SNP_Base") == pl.col("base_code")).alias("Match")
        )

        het_df = df.filter(~pl.col("base_code").is_in(fixed_codes)).with_columns(
            degenerate_match_expr()
        )

        ploidy_df = df.filter(pl.col("base_code") < 0)

        fixed_n  = fixed_df.shape[0]
        fixed_m  = fixed_df.filter(pl.col("Match")).shape[0]

        het_n    = het_df.shape[0]
        het_m    = het_df.filter(pl.col("Match")).shape[0]

        ploidy_n = ploidy_df.shape[0]

        total_n  = fixed_n + het_n + ploidy_n

        fixed_prop = (fixed_m / fixed_n) if fixed_n > 0 else 0.0
        het_prop   = (het_m / het_n)     if het_n > 0 else 0.0

        sample_rows.append([
            sample_id, gr,
            total_n,
            fixed_n, fixed_m, fixed_prop,
            het_n, het_m, het_prop,
            ploidy_n,
        ])

    result = pl.DataFrame(
        sample_rows,
        schema=[
            "Sample_ID",
            "Focal_Group",
            "Covered",
            "Fixed",
            "Fixed_Match",
            "Fixed_Match_Prop",
            "Heterozygous",
            "Het_Match",
            "Het_Match_Prop",
            "Ploidy_Fail",
        ],
        orient="row"
    )

    return result


def getCalledBases(path):

    if not os.path.exists(path):
        raise ValueError(f"Called bases path does not exist: {path}")

    if os.path.isfile(path) and path.endswith("_Called.parquet"):
        return [os.path.abspath(path)]

    if os.path.isfile(path) and not path.endswith("_Called.parquet"):
        with open(path) as f:
            files = [line.strip() for line in f if line.strip()]

        if not files:
            raise ValueError(f"List file is empty: {path}")

        bad = [f for f in files if not f.endswith("_Called.parquet")]

        if bad:
            raise ValueError(
                "List file contains entries not ending in '_Called.parquet':\n" +
                "\n".join(bad)
            )

        missing = [f for f in files if not os.path.exists(f)]
        if missing:
            raise ValueError(
                "Some paths in the list file do not exist:\n" +
                "\n".join(missing)
            )

        return [os.path.abspath(file) for file in files]

    raise ValueError(
        f"`--called_bases` must be either:\n"
        f"  • A *_Called.parquet file, OR\n"
        f"  • A text file listing *_Called.parquet files\n"
        f"Given: {path}"
    )

args = parser.parse_args()

called_base_paths = getCalledBases(args.called_bases)
snp_df = pl.read_parquet(args.snp_parquet)

called_results = []

for called_path in called_base_paths:

    sample_id = os.path.basename(called_path).replace("_Called.parquet","")
    called_results.append(snpClassifier(snp_df,sample_id,called_path))

called_class_results = pl.concat(called_results)

if not args.output_file.endswith(".tsv"):
    output_file = os.path.abspath(f"{args.output_file}.tsv")
else:
    output_file = os.path.abspath(f"{args.output_file}")
    
called_class_results.to_pandas().to_csv(output_file,sep="\t")