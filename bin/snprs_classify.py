#!/usr/bin/env python3
import pandas as pd
import os
import sys
from natsort import natsorted
from collections import Counter
import argparse
import polars as pl
import pyarrow.parquet as pq
import shutil
import json

def parse_args():
    parser = argparse.ArgumentParser(description="Classify SNPRS called base info based on a diagnostic SNP dataset")
    
    # Data args
    parser.add_argument("--out", dest="output_directory", type=str, required=True,help="Path to save output")
    parser.add_argument("--sample_id", dest="sample_id", type=str, required=True,help="Sample ID")
    parser.add_argument("--called_base", dest="called_base_file", type=str, required=True,help="Path to _Called_Bases.parquet")
    parser.add_argument("--group_file", dest="group_file", type=str, required=True,help="Path to group file used to generate SNPs")
    parser.add_argument("--comparison_file", dest="comparison_file", type=str, required=True,help="Path to _Comparisons.csv")
    parser.add_argument("--snp_parquet", dest="snp_parquet", type=str, required=True,help="Path to _SNPs.parquet")
    parser.add_argument("--half",dest="require_half",action="store_true",help="Require half of ingroup/outgroup to call a SNP")
    return parser.parse_args()


pl.Config.set_tbl_rows(-1)
pl.Config.set_tbl_cols(-1)

# region 00: Parse args
args = parse_args()

output_directory = os.path.abspath(args.output_directory)
sample_id = str(args.sample_id)
called_base_file = os.path.abspath(args.called_base_file)
snp_parquet_file = os.path.abspath(args.snp_parquet)
group_file = os.path.abspath(args.group_file)
comparison_file = os.path.abspath(args.comparison_file)

if not os.path.exists(called_base_file):
    sys.exit(f"${called_base_file} does not exist...")
elif(not os.path.exists(snp_parquet_file)):
    sys.exit(f"${snp_parquet_file} does not exist...")
elif(not os.path.exists(group_file)):
    sys.exit(f"${group_file} does not exist...")
elif(not os.path.exists(comparison_file)):
    sys.exit(f"${comparison_file} does not exist...")
elif not os.path.exists(output_directory):
    parent = os.path.dirname(output_directory)
    if os.path.exists(parent):
        os.makedirs(output_directory)
    else:
        sys.exit(f"Cannot create {output_directory}...")

# Group data
group_df = pd.read_csv(group_file, usecols=["Clade_ID", "Clade_Tips", "Clade_Type"]).copy()
group_df["Tip_Count"] = group_df["Clade_Tips"].str.count(";") + 1
clade_tip_dict = dict(zip(group_df["Clade_ID"], group_df["Tip_Count"]))

terminal_df  = group_df.query('Clade_Type == "Terminal"').copy()
internal_df  = group_df.query('Clade_Type == "Internal"').copy()

terminal_ids = list(terminal_df['Clade_ID'])
internal_ids = list(internal_df['Clade_ID'])

# Comparison data
comparison_df = pd.read_csv(comparison_file, usecols=["Ingroup", "Comparison", "SNP_Count"]).query("SNP_Count > 0").copy()
comparison_pl = pl.from_pandas(comparison_df)

# SNP data
lazy_snp_data = pl.scan_parquet(snp_parquet_file).rename({"Outgroup":"Comparison"})

if args.require_half:

    in_half_map = pl.LazyFrame({
        "Ingroup": list(clade_tip_dict.keys()),
        "In_Tip_Count": list(clade_tip_dict.values())
    })

    lazy_snp_data = (
        lazy_snp_data
        .join(in_half_map, on="Ingroup", how="left")
        .with_columns([
            (pl.col("In_Tip_Count") // 2).alias("In_Half_Tips")
        ])
        .filter(
            (pl.col("In_Count") >= pl.col("In_Half_Tips"))
            )
        .drop(["In_Tip_Count","In_Half_Tips"])
    )

lazy_snp_count = (
    lazy_snp_data
    .group_by(["Ingroup", "Comparison"])
    .agg(pl.len().alias("SNP_Count"))
).collect()

# endregion

# region 01: Get covered sites

lazy_called_bases = pl.scan_parquet(called_base_file).filter(pl.col('type').is_in([0,1,3,4]))

joined_called_bases = lazy_snp_data.join(lazy_called_bases,on=['contig_index','contig_position'],how="left")

type_map = {
    "0": "Fixed_Base",
    "1": "Fixed_Deletion",
    "3": "Het_Base",
    "4": "Het_Deletion",
    "null": "Missing"
}

type_counts = (
    joined_called_bases
    .unique(["contig_index", "contig_position"])
    .select(pl.col("type").cast(pl.Utf8))
    .with_columns(
        pl.col("type").fill_null("null").replace_strict(type_map).alias("Site_Type")
    )
    .select("Site_Type")
    .group_by("Site_Type")
    .agg(pl.len().alias("Type_Count"))
)

lazy_covered_snps = (
    joined_called_bases
    .filter(pl.col("type").is_not_null())
)

# endregion

# region 02: Simple rows

fixed_bases = {'A','C','T','G','-','?'}

simple_rows = (
    lazy_covered_snps
    .filter(
        pl.col("final_base").is_in(fixed_bases) &
        pl.col("In_Base").is_in(fixed_bases)
    )
    .with_columns(
        (pl.col("final_base") == pl.col("In_Base"))
        .cast(pl.Int8)
        .alias("Match")
    )
)
simple_row_count = simple_rows.select(pl.len()).collect().item()

simple_match_counts = pl.DataFrame()

if simple_row_count > 0:

    simple_match_counts = (
        simple_rows
        .select("Ingroup", "Comparison", "Match")
        .group_by("Ingroup", "Comparison", "Match")
        .agg(pl.len().alias("Match_Count"))
        .collect()
        .pivot(
            index=["Ingroup", "Comparison"],
            on="Match",
            values="Match_Count"
        )
        .fill_null(0)
        .rename({"0": "No_Match", "1": "Match"})
        .with_columns([
            (pl.col("Match") + pl.col("No_Match")).alias("Covered"),
            pl.when((pl.col("Match") + pl.col("No_Match")) > 0)
            .then(pl.col("Match") / (pl.col("Match") + pl.col("No_Match")))
            .otherwise(0)
            .alias("Ratio")
        ])
    )

# endregion

# region 02: Degenerate sample rows

degenerate_map = {
"A": ["A"],
"C": ["C"],
"G": ["G"],
"T": ["T"],
"-": ["-"],
"R": ["A","G"],
"Y": ["C","T"],
"S": ["G","C"],
"W": ["A","T"],
"K": ["G","T"],
"M": ["A","C"],
"B": ["C","G","T"],
"D": ["A","G","T"],
"H": ["A","C","T"],
"V": ["A","C","G"],

"a": ["A","-"],
"c": ["C","-"],
"g": ["G","-"],
"t": ["T","-"],
"r": ["A","G","-"],
"y": ["C","T","-"],
"s": ["G","C","-"],
"w": ["A","T","-"],
"k": ["G","T","-"],
"m": ["A","C","-"],
"b": ["C","G","T","-"],
"d": ["A","G","T","-"],
"h": ["A","C","T","-"],
"v": ["A","C","G","-"],
}

degenerate_sample_rows = (
    lazy_covered_snps
    .filter(
        pl.col("In_Base").is_in(fixed_bases)
        & (~pl.col("final_base").is_in(fixed_bases))
    )
    .with_columns([
        pl.col("final_base")
        .map_elements(lambda v: "".join(sorted(set(degenerate_map.get(v, [])))), return_dtype=pl.String)
        .alias("sample_bases")
    ])
    .with_columns(
        pl.struct(["In_Base", "sample_bases"])
        .map_elements(lambda x: 1 if x["In_Base"] in x["sample_bases"] else 0,return_dtype=pl.Int8)
        .alias("Match")
        .cast(pl.Int8)
    )
    .drop('sample_bases')
)

degenerate_row_count = degenerate_sample_rows.select(pl.len()).collect().item()

degenerate_match_counts = pl.DataFrame()

if degenerate_row_count > 0:

    degenerate_match_counts = (
        degenerate_sample_rows
        .select("Ingroup", "Comparison", "Match")
        .group_by("Ingroup", "Comparison", "Match")
        .agg(pl.len().alias("Match_Count"))
        .collect()
        .pivot(
            index=["Ingroup", "Comparison"],
            on="Match",
            values="Match_Count"
        )
        .fill_null(0)
        .rename({"0": "No_Match", "1": "Match"})
        .with_columns([
            (pl.col("Match") + pl.col("No_Match")).alias("Covered"),
            pl.when((pl.col("Match") + pl.col("No_Match")) > 0)
            .then(pl.col("Match") / (pl.col("Match") + pl.col("No_Match")))
            .otherwise(0)
            .alias("Ratio")
        ])
    )

# region 03: Process Ingroup comparisons

ingroup_comparisons = (
    pl.concat([
        simple_rows.with_columns(pl.lit("Simple").alias("Source")),
        degenerate_sample_rows.with_columns(pl.lit("Degenerate").alias("Source")),
    ])
    .filter(pl.col("Comparison") == "Ingroup")
)

ingroup_count_summary = (
    ingroup_comparisons
    .unique(['contig_index','contig_position'])
    .group_by("Ingroup")
    .agg([
        pl.col("In_Count").quantile(0.25).alias("In_Q25"),
        pl.col("In_Count").quantile(0.50).alias("In_Q50"),
        pl.col("In_Count").quantile(0.75).alias("In_Q75"),
        pl.col("Out_Count").quantile(0.25).alias("Out_Q25"),
        pl.col("Out_Count").quantile(0.50).alias("Out_Q50"),
        pl.col("Out_Count").quantile(0.75).alias("Out_Q75"),
        pl.col("Out_Groups").quantile(0.25).alias("Outg_Q25"),
        pl.col("Out_Groups").quantile(0.50).alias("Outg_Q50"),
        pl.col("Out_Groups").quantile(0.75).alias("Outg_Q75")
    ])
    .select(['Ingroup','In_Q25','In_Q50','In_Q75','Out_Q25','Out_Q50','Out_Q75','Outg_Q25','Outg_Q50','Outg_Q75'])
)


ingroup_comparisons = (
    ingroup_comparisons
    .join(ingroup_count_summary, on="Ingroup")
    .with_columns([
        pl.when(pl.col("In_Count") < pl.col("In_Q25"))
        .then(pl.lit("Low"))
        .when(pl.col("In_Count") >= pl.col("In_Q75"))
        .then(pl.lit("High"))
        .otherwise(pl.lit("Moderate"))
        .alias("Ingroup_Count_Confidence"),
        pl.when(pl.col("Out_Count") < pl.col("Out_Q25"))
        .then(pl.lit("Low"))
        .when(pl.col("Out_Count") >= pl.col("Out_Q75"))
        .then(pl.lit("High"))
        .otherwise(pl.lit("Moderate"))
        .alias("Outgroup_Count_Confidence"),
        pl.when(pl.col("Out_Groups") < pl.col("Outg_Q25"))
        .then(pl.lit("Low"))
        .when(pl.col("Out_Groups") >= pl.col("Outg_Q75"))
        .then(pl.lit("High"))
        .otherwise(pl.lit("Moderate"))
        .alias("Outgroup_Number_Confidence")]
    )
    .select([
    'contig_index','contig_position','type',
    'Ingroup','Comparison',
    'In_Count','Ingroup_Count_Confidence',
    'Out_Count','Outgroup_Count_Confidence',
    'Out_Groups','Outgroup_Number_Confidence',
    'Match','Source'
    ])
    .filter((pl.col("Ingroup_Count_Confidence")!="Low") & (pl.col("Outgroup_Count_Confidence")!="Low")  & (pl.col("Outgroup_Count_Confidence")!="Low"))
)

ingroup_match_counts = (
    ingroup_comparisons
    .select("Ingroup", "Comparison", "Match")
    .group_by("Ingroup", "Comparison", "Match")
    .agg(pl.len().alias("Match_Count"))
    .collect()
    .pivot(
        index=["Ingroup", "Comparison"],
        on="Match",
        values="Match_Count"
    )
    .fill_null(0)
    .rename({"0": "No_Match", "1": "Match"})
    .with_columns([
        (pl.col("Match") + pl.col("No_Match")).alias("Covered"),
        pl.when((pl.col("Match") + pl.col("No_Match")) > 0)
        .then(pl.col("Match") / (pl.col("Match") + pl.col("No_Match")))
        .otherwise(0)
        .alias("Ratio")
    ])
    .sort('Ratio')
)
print(ingroup_match_counts.filter((pl.col("Ingroup").is_in(terminal_ids) & (pl.col("Ratio") > 0.80))))
print(ingroup_match_counts.filter((pl.col("Ingroup").is_in(internal_ids) & (pl.col("Ratio") > 0.80))))

"""

    .with_columns(
        pl.when(pl.col("Out_Count") < pl.col("Out_Q25"))
          .then("Low")
          .when(pl.col("Out_Count") > pl.col("Out_Q75"))
          .then("High")
          .otherwise("Moderate")
          .alias("Confidence_Level")
    )
)



ingroup_df   = group_df.query('Clade_Type == "Ingroup"').copy()
outgroup_df  = group_df.query('Clade_Type == "Outgroup"').copy()

for terminal in terminal_ids:

    internal_terminal_df = pl.DataFrame()
    terminal_ingroup_df = pl.DataFrame()
    terminal_internal_df = pl.DataFrame()


    terminal_count_df = (
        simple_match_counts
        .filter(pl.col("Ingroup") == terminal)
    )

    terminal_ingroup_df = (
        terminal_count_df
        .filter(pl.col("Comparison") == "Ingroup")
        .select(["Ingroup", "Comparison","Covered","No_Match","Match", "Ratio"])
        .join(lazy_snp_count,on=["Ingroup", "Comparison"],how="left")
        .with_columns((pl.col("Covered") / pl.col("SNP_Count")).alias("Percent_Covered"))
        .select("Ingroup","Comparison","SNP_Count","Covered","Percent_Covered",'Match','Ratio')
        .sort('Ratio')
    )

    comparison_list = comparison_df[(comparison_df['Ingroup'] == terminal) & (comparison_df['Comparison'] != "Ingroup")]['Comparison'].tolist()
    
    if len(comparison_list) > 0:

        terminal_internal_df = (
            terminal_count_df
            .filter(pl.col("Comparison").is_in(comparison_list))
            .select(["Ingroup", "Comparison","Covered","No_Match","Match", "Ratio"])
            .join(lazy_snp_count,on=["Ingroup", "Comparison"],how="left")
            .with_columns((pl.col("Covered") / pl.col("SNP_Count")).alias("Percent_Covered"))
            .select("Ingroup","Comparison","SNP_Count","Covered","Percent_Covered",'Match','Ratio')
            .sort('Ratio')
        )

        internal_terminal_df = (
            simple_match_counts
            .filter(pl.col("Ingroup").is_in(comparison_list))
            .select(["Ingroup", "Comparison","Covered","No_Match","Match", "Ratio"])
            .join(lazy_snp_count,on=["Ingroup", "Comparison"],how="left")
            .with_columns((pl.col("Covered") / pl.col("SNP_Count")).alias("Percent_Covered"))
            .select("Ingroup","Comparison","SNP_Count","Covered","Percent_Covered",'Match','Ratio')
            .sort('Ratio')
        )
        
    full_df = pl.concat([terminal_ingroup_df,terminal_internal_df,internal_terminal_df])
    
    terminal_full_df = full_df.filter(pl.col("Ingroup") == terminal)

    print(terminal)
    print(full_df.filter(pl.col('Ratio')>0.9))



for term_id in terminal_ids:
    term_ingroup = simple_match_counts.filter((pl.col('Ingroup') == term_id)).sort('Outgroup_Ratio','Ratio')
    print(term_ingroup)

# region XX: Degenerate inbase rows?
# endregion
"""