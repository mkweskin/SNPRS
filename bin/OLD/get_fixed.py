import polars as pl
import os
import pyarrow.parquet as pq
import json 
import argparse
import sys
import math
from datetime import datetime

def process_fixed(fixed_id,scaffold_file,base_file,code_file,site_file,gap_string,group_samples,missing_allowed,output_directory,joined_directory,ref_fasta):
    
    output_parquet = os.path.join(output_directory,f"{fixed_id}.parquet")

    must_have_data = [sample for sample in group_samples if sample not in missing_allowed]

    scaffold_df = pl.read_parquet(scaffold_file).lazy()
    code_df = pl.read_parquet(code_file).lazy().select(group_samples)
    base_df = pl.read_parquet(base_file).lazy().select(group_samples)

    combined_code_df = pl.concat([scaffold_df, code_df], how="horizontal")
    combined_base_df = pl.concat([scaffold_df, base_df], how="horizontal")

    if gap_string == "1":
        valid_bases = {"A", "C", "T", "G", "-"}
    else:
        valid_bases = {"A", "C", "T", "G"} 

    # Get sites that are non-singleton if 2+, or fixed ACTG- if 1 sample
    base_filters = []
    
    for sample in must_have_data:
        base_filters.append(pl.col(sample).is_in(valid_bases))
    for sample in missing_allowed:
        base_filters.append(pl.col(sample).is_in(valid_bases.union({"N"})))
    
    if len(group_samples) > 1:

        code_filters = []
        for sample in must_have_data:
            code_filters.append(pl.col(sample) == 1)
        for sample in missing_allowed:
            code_filters.append(pl.col(sample).is_in([0, 1]))
        
        filtered_code_df = combined_code_df.filter(pl.all_horizontal(*code_filters)).select(['contig_index','contig_position'])
        filtered_base_df = filtered_code_df.join(combined_base_df, on=['contig_index','contig_position'], how="left").filter(pl.all_horizontal(*base_filters)).collect()
        
        final_base_df = (
            filtered_base_df
            .filter(
                pl.struct(group_samples)
                .map_elements(lambda row: len({v for v in row.values() if v != "N"}), return_dtype=pl.Int8) == 1
            )
            .with_columns(
                pl.struct(group_samples)
                .map_elements(lambda row: next(v for v in row.values() if v != "N"), return_dtype=pl.String)
                .alias(fixed_id)
            )
            .select(['contig_index', 'contig_position', fixed_id])
        )

    else:
        sample_id = group_samples[0]
        final_base_df = (
            combined_base_df
            .filter(pl.col(sample_id).is_in(valid_bases))
            .select(['contig_index','contig_position',sample_id])
            .rename({sample_id: fixed_id})
            .collect()
        )
    
    site_count = final_base_df.height

    metadata = {
        "Joined_Directory":joined_directory,
        "Grouped_Sample":",".join(group_samples),
        "Ref_FASTA":ref_fasta,
        "Sample_Count":str(len(group_samples)),
        "Gaps_Included":"True" if gap_string == 1 else "False",
        "Allowed_Missing":",".join(missing_allowed),
        "Final_Site_Count":str(site_count)
    }

    final_metadata = {k.encode(): v.encode() for k, v in metadata.items()}
    final_base_df = final_base_df.to_arrow().replace_schema_metadata(final_metadata)
    pq.write_table(final_base_df, output_parquet, compression="snappy")
        
if __name__ == "__main__":

    fixed_id,scaffold_file,base_file,code_file,site_file,gap_string,group_file,missing_file,output_directory,joined_directory,ref_fasta = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10],sys.argv[11]

    # Get group and missing info again
    group_samples = []
    with open(group_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            group_samples.append(line)

    missing_allowed = []
    if str(missing_file) != "None":
        with open(missing_file, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                missing_allowed.append(line)
    
    process_fixed(fixed_id,scaffold_file,base_file,code_file,site_file,gap_string,group_samples,missing_allowed,output_directory,joined_directory,ref_fasta)
    