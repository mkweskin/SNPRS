#!/usr/bin/env python3

import polars as pl
import os
import argparse


parser = argparse.ArgumentParser(description='ID called bases from a SNP parquet')
parser.add_argument('-r','--raw_parquets',dest="raw_parquets", type=str,required=True, help='Path to a SNPRS raw parquet file or a file with multiple paths.')
parser.add_argument('-s','--snp',dest="snp_parquet", type=str,required=True, help='Path to SNP parquet file')
parser.add_argument('-n','--name',dest="output_file", default=None,type=str, help='Name for output file')

def rawSNPClassifer(snp_df, sample_id, raw_parquet_path):

    convert_dict = {1:"A",2:"C",3:"G",4:"T",16:"-"}

    snp_groups = snp_df["SNP_Group"].unique().to_list()
    snp_coords = snp_df.select(['contig_index','contig_position']).unique()
    
    sample_parquet = pl.read_parquet(raw_parquet_path).join(snp_coords,on=['contig_index','contig_position'],how="inner")
    
    sample_rows = []

    for sp in snp_groups:
        
        sp_snp_df = snp_df.filter((pl.col("SNP_Group")==sp))

        match_df = (
            sp_snp_df.join(sample_parquet, on=['contig_index','contig_position'], how="left")
            .select(['contig_index','contig_position','SNP_Base','base','frequency','depth'])
            .filter(~pl.col('depth').is_null())
            .with_columns([
                pl.col('SNP_Base').replace_strict(convert_dict).alias('SNP_Base')
            ])
            .with_columns([
                (pl.col('SNP_Base') == pl.col("base")).alias("Match")
            ])
            .with_columns([
                (pl.col('frequency') * pl.col('depth')).alias("allele_depth")
            ])
        )

        distinct_count = (
            match_df
            .group_by(["contig_index", "contig_position"])
            .len()
            .height
        )

        match_summary_df = (
            match_df
            .group_by(['Match'])
            .agg([
                pl.col("allele_depth").sum().alias("Sum_Allele_Depth")
            ])
        )

        count_match_df = match_summary_df.filter(pl.col("Match") == True)
        count_nonmatch_df = match_summary_df.filter(pl.col("Match") == False)

        if count_match_df.height == 0:
            sum_match = 0
        else:
            sum_match = count_match_df.select("Sum_Allele_Depth").item()

        if count_nonmatch_df.height == 0:
            sum_nonmatch = 0
        else:
            sum_nonmatch = count_nonmatch_df.select("Sum_Allele_Depth").item()

        sample_rows.append([sample_id,sp,distinct_count,sum_match, sum_nonmatch])

    sample_df = pl.DataFrame(
        sample_rows,
        schema=[
            "Sample_ID",
            "Reference_Species",
            "SNP_Count",
            "Match",
            "Non_Match"
        ],
        orient="row"
    ).with_columns([
    pl.col("SNP_Count").cast(pl.Int32),
    pl.col("Match").cast(pl.Int32),
    pl.col("Non_Match").cast(pl.Int32),
])

    return sample_df

def getRawParquets(path):

    if not os.path.exists(path):
        raise ValueError(f"Raw parquet path does not exist: {path}")

    if os.path.isfile(path) and path.endswith("_Raw.parquet"):
        return [os.path.abspath(path)]

    if os.path.isfile(path) and not path.endswith("_Raw.parquet"):
        with open(path) as f:
            files = [line.strip() for line in f if line.strip()]

        if not files:
            raise ValueError(f"List file is empty: {path}")

        bad = [f for f in files if not f.endswith("_Raw.parquet")]

        if bad:
            raise ValueError(
                "List file contains entries not ending in '_Raw.parquet':\n" +
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
        f"`--raw_parquets` must be either:\n"
        f"  • A *_Raw.parquet file, OR\n"
        f"  • A text file listing *_Raw.parquet files\n"
        f"Given: {path}"
    )

args = parser.parse_args()

raw_parquet_paths = getRawParquets(args.raw_parquets)
snp_df = pl.read_parquet(args.snp_parquet)

called_results = []

for raw_path in raw_parquet_paths:

    sample_id = os.path.basename(raw_path).replace("_Raw.parquet","")
    print(sample_id)
    called_results.append(rawSNPClassifer(snp_df,sample_id,raw_path))

called_class_results = pl.concat(called_results)

if not args.output_file.endswith(".tsv"):
    output_file = os.path.abspath(f"{args.output_file}.tsv")
else:
    output_file = os.path.abspath(f"{args.output_file}")
    
called_class_results.to_pandas().to_csv(output_file,sep="\t")