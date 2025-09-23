import polars as pl
import json
import sys
import os

def filter_missing(site_file,pass_site_type_file,max_missing,temp_directory):
    
    valid_missing_file = os.path.join(temp_directory, "Valid_Missing.json")

    with open(pass_site_type_file, "r", encoding="utf-8") as f:
        pass_site_type_rows = set(json.load(f))

    sites = pl.scan_parquet(site_file).with_row_index("row_nr")
    mask = pl.col("row_nr").is_in(list(pass_site_type_rows))
    
    missing_mask = pl.lit(False)
    missing_mask |= pl.col("missing") <= int(max_missing)

    # Combine gross exclusion mask with type mask
    final_mask = mask & missing_mask

    surviving_rows = sites.filter(final_mask).collect().select("row_nr")
    surviving_row_numbers = surviving_rows["row_nr"].to_list()
    
    if len(surviving_row_numbers) == 0:
        sys.exit("No data remains after filtering for requested site types")
        
    with open(valid_missing_file, "w", encoding="utf-8") as f:
        json.dump(surviving_row_numbers, f)

    print(",".join([valid_missing_file, str(len(surviving_row_numbers))]))

if __name__ == "__main__":

    site_file,pass_site_type_file,max_missing,temp_directory = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]
    filter_missing(site_file,pass_site_type_file,max_missing,temp_directory)
