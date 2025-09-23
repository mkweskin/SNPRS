import sys
import polars as pl
import os
import json
from collections import defaultdict, Counter
from natsort import natsorted
from concurrent.futures import ProcessPoolExecutor
import math

def classify_row(row):

    valid_bases = {'A', 'C', 'G', 'T', '-'}
    fixed_counts = Counter(base for base in row if base in valid_bases)
    singleton_bases = {b for b, c in fixed_counts.items() if c == 1}

    counts = {
        "missing": 0,
        "heterozygous": 0,
        "singleton": 0,
        "nonsingleton": 0,
        "invalid":0,
        "gap":0
    }
    
    nonsingleton_alleles = set()
    codes = []

    for base in row:
        if base is None:
            codes.append(0)
            counts["missing"] += 1
            
        elif base not in valid_bases:
            if base == "N":
                codes.append(4)
                counts["invalid"] += 1
            else:
                codes.append(3)
                counts["heterozygous"] += 1
        
        elif base in singleton_bases:
            codes.append(2)
            counts["singleton"] += 1
            if base == "-":
                counts["gap"] += 1
        
        else:
            codes.append(1)
            counts["nonsingleton"] += 1
            nonsingleton_alleles.add(base)
            if base == "-":
                counts["gap"] += 1
    
    missing_count = counts["missing"] + counts['invalid']
    singleton_count = counts["singleton"]
    allele_count = len(nonsingleton_alleles)
    counts["nonsingleton_alleles"] = allele_count

    # Site codes
    # 0: All singletons: 0 nonsingleton alleles
    # 1: Pure fixed: 1 allele, 0 singletons
    # 2: Pure biallelic: 2 alleles, 0 singletons 
    # 3: Pure triallelic: 3 alleles, 0 singletons 
    # 4: Pure quadallelic: 4 alleles, 0 singletons 
    # 5: Pure pentallelic: 5 alleles, 0 singletons 
    # 6: Fixed w/singletons: 1 allele, 1+ singletons
    # 7: Biallelic w/singletons: 2 alleles, 1+ singletons 
    # 8: Triallelic w/singletons: 3 alleles, 1+ singletons 
    # 9: Quadallelic w/singletons: 4 alleles, 1+ singletons 
    # 10: Pentallelic w/singletons: 5 alleles, 1+ singletons (???)
    
    def get_missing_dict(missing_count, singleton_count, allele_count):
        if singleton_count == 0:
            return {missing_count:allele_count}
        elif allele_count == 0:
            return {missing_count:0}
        else:
            return {missing_count:5 + allele_count}
        
    missing_dict = get_missing_dict(missing_count,singleton_count,allele_count)
    return counts, codes, missing_dict

def populate_chunk(sample_chunk_paths,i,temp_directory,sample_ids):
    
    dfs = [pl.scan_parquet(path).select("final_base").rename({"final_base": sample_id}) for path, sample_id in zip(sample_chunk_paths, sample_ids)]
    df = pl.concat(dfs, how="horizontal").collect()
    results = [classify_row(row) for row in df.iter_rows()]

    counts, codes, missing = zip(*results)
    
    site_df = pl.DataFrame(counts,orient="row")
    code_df = pl.DataFrame(codes, schema=sample_ids,orient="row")
    
    site_file = os.path.join(temp_directory, f"Sites_Chunk_{i}.parquet")
    code_file = os.path.join(temp_directory, f"Codes_Chunk_{i}.parquet")

    site_df.write_parquet(site_file, compression="snappy")
    code_df.write_parquet(code_file, compression="snappy")
    
    grouped_missing = defaultdict(Counter)

    for d in missing:
        for missing_count, code in d.items():
            grouped_missing[missing_count][code] += 1

    grouped_missing = {k: dict(v) for k, v in grouped_missing.items()}    
    
    output = {
    "site_file": site_file,
    "code_file": code_file,
    "missing_data":grouped_missing,
    "chunk_id":str(i)
    }
    
    print(json.dumps(output))

if __name__ == "__main__":
    chunk_file,temp_directory = sys.argv[1], sys.argv[2]
    
    with open(chunk_file, "r") as f:
        chunk_list = json.load(f)
    
    sample_id_paths = [(os.path.basename(path).split("_Called")[0], path) for path in chunk_list]
    sorted_path_sample_pairs = natsorted(sample_id_paths, key=lambda x: x[0])
    sample_ids, chunk_list_sorted = zip(*sorted_path_sample_pairs)

    chunk_number = int(os.path.splitext(os.path.basename(chunk_list_sorted[0]).split("_Chunk_")[-1])[0])

    populate_chunk(chunk_list_sorted,chunk_number,temp_directory,sample_ids)
