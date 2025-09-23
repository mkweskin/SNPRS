import polars as pl
import json
import sys
import os

def filter_valid_site_types(site_file,pass_gross_exclusion_file,fixed_string,bi_string,tri_string,quad_string,pent_string,sing_string,temp_directory):
    
    valid_sites_file = os.path.join(temp_directory, "Valid_Site_Types.json")

    with open(pass_gross_exclusion_file, "r", encoding="utf-8") as f:
        gross_exclusion_rows = set(json.load(f))

    sites = pl.scan_parquet(site_file).with_row_index("row_nr")
    mask = pl.col("row_nr").is_in(list(gross_exclusion_rows))
    type_mask = pl.lit(False)

    if bi_string == "1":
        type_mask |= pl.col("nonsingleton_alleles") == 2
    if tri_string == "1":
        type_mask |= pl.col("nonsingleton_alleles") == 3
    if quad_string == "1":
        type_mask |= pl.col("nonsingleton_alleles") == 4
    if pent_string == "1":
        type_mask |= pl.col("nonsingleton_alleles") == 5
    if sing_string == "1":
        type_mask |= (pl.col("singleton") >= 1) & (pl.col("nonsingleton_alleles") == 1)
    if fixed_string == "1":
        type_mask |= (pl.col("singleton") == 0) & (pl.col("nonsingleton_alleles") == 1)

    # Combine gross exclusion mask with type mask
    final_mask = mask & type_mask

    surviving_rows = sites.filter(final_mask).collect().select("row_nr")
    surviving_row_numbers = surviving_rows["row_nr"].to_list()
    
    if len(surviving_row_numbers) == 0:
        sys.exit("No data remains after filtering for requested site types")
        
    with open(valid_sites_file, "w", encoding="utf-8") as f:
        json.dump(surviving_row_numbers, f)

    print(",".join([valid_sites_file, str(len(surviving_row_numbers))]))

if __name__ == "__main__":

    site_file,pass_gross_exclusion_file,fixed_string,bi_string,tri_string,quad_string,pent_string,sing_string,temp_directory = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9]
    filter_valid_site_types(site_file,pass_gross_exclusion_file,fixed_string,bi_string,tri_string,quad_string,pent_string,sing_string,temp_directory)
