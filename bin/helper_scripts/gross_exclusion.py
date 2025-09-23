import polars as pl
import json
import sys
import os

def gross_exclusion(site_file, gap_string, het_string, invalid_string, sing_string, temp_directory):
    gross_exclusion_file = os.path.join(temp_directory, "Gross_Exclusion.json")

    sites = pl.scan_parquet(site_file).with_row_index("row_nr")
    mask = pl.lit(True)

    if het_string == "0":
        mask &= pl.col("heterozygous") == 0

    if sing_string == "0":
        mask &= pl.col("singleton") == 0

    if gap_string == "0":
        mask &= pl.col("gap") == 0

    if invalid_string == "0":
        mask &= pl.col("invalid") == 0
    
    surviving_rows = sites.filter(mask).collect().select("row_nr")
    surviving_row_numbers = surviving_rows["row_nr"].to_list()

    if len(surviving_row_numbers) == 0:
        sys.exit("No data remains after removing het/gap/invalid/singetons")
        
    with open(gross_exclusion_file, "w", encoding="utf-8") as f:
        json.dump(surviving_row_numbers, f)

    print(",".join([gross_exclusion_file, str(len(surviving_row_numbers))]))

if __name__ == "__main__":

    site_file,gap_string,het_string,invalid_string,sing_string,temp_directory = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6]
    gross_exclusion(site_file,gap_string,het_string,invalid_string,sing_string,temp_directory)
