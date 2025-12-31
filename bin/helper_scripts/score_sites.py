import sys
import polars as pl
import os
import json
from collections import defaultdict, Counter
from natsort import natsorted
from concurrent.futures import ProcessPoolExecutor,as_completed
import math
import argparse
import pyarrow.parquet as pq
import glob
import multiprocessing
import time

fixed_codes = {1, 2, 3, 4, 16}
het_codes = set(range(5, 16)) | set(range(17, 32))
cpu = multiprocessing.cpu_count()

def parse_args():
    parser = argparse.ArgumentParser(description="Create scaffold parquet from called base files")
    
    parser.add_argument("--base", dest="base_file", type=str, required=True,help="Path to base parquet")
    parser.add_argument("--scaffold", dest="scaffold_file", type=str, required=True,help="Path to scaffold parquet")
    parser.add_argument("--out_dir", dest="output_directory", type=str, required=True,help="Path to output Parquet files")
    parser.add_argument("--join_id", dest="join_id", type=str, required=True,help="Prefix for output files")
    parser.add_argument("--mem_mode", dest="mem_mode", action = "store_true", help="Split chunks up 10X smaller than default")

    return parser.parse_args()

def process_chunk(i, start, stop, base_file, temp_directory):
    t0 = time.time()
    #print(f"[{i}] Starting process_chunk: {time.ctime()}")

    # 1️⃣ Fetch sample IDs
    sample_ids = pl.scan_parquet(base_file).collect_schema().names()

    # 2️⃣ Set up output files
    code_chunk_file = os.path.join(temp_directory, f"Code_Chunk_{i}.parquet")
    site_chunk_file = os.path.join(temp_directory, f"Site_Chunk_{i}.parquet")

    length = stop - start

    # 3️⃣ Collect DataFrame slice
    t1 = time.time()
    df_chunk = pl.scan_parquet(base_file).slice(start, length).collect()

    # 4️⃣ Convert rows
    rows = list(df_chunk.iter_rows())

    # 5️⃣ Classify rows
    results = [classify_row(list(row)) for row in rows]

    # 6️⃣ Unpack results
    counts_list, codes_list, missing_dicts = zip(*results)

    # 7️⃣ Write Parquet files
    pl.DataFrame(counts_list, orient="row").write_parquet(site_chunk_file, compression="snappy")
    pl.DataFrame(codes_list, schema=sample_ids, orient="row").write_parquet(code_chunk_file, compression="snappy")

    # 8️⃣ Aggregate missing counts
    grouped_missing = defaultdict(Counter)
    for d in missing_dicts:
        for missing_count, code in d.items():
            grouped_missing[missing_count][code] += 1
    grouped_missing = {k: dict(v) for k, v in grouped_missing.items()}

    total_time = time.time() - t0
    #print(f"[{i}] Finished process_chunk in {total_time:.2f}s ({time.ctime()})\n")

    return code_chunk_file, site_chunk_file, grouped_missing

def classify_row(row):

    counts = Counter(row) 
    
    missing_count = counts.get(0, 0)
    ploidy_fail_count = sum(c for k, c in counts.items() if k < 0)
    het_count = sum(c for b, c in counts.items() if b in het_codes)

    fixed_counts = {b: c for b, c in counts.items() if b in fixed_codes}

    if fixed_counts:
        
        nonsingleton_bases = {b:c for b, c in fixed_counts.items() if c > 1}

        if nonsingleton_bases:
            nonsingleton_sample_count = sum(c for b, c in nonsingleton_bases.items())
            nonsingleton_alleles = {b for b,c in nonsingleton_bases.items()}
            nonsingleton_allele_count = len(nonsingleton_alleles)
            row = [(32 + b) if b in nonsingleton_alleles else b for b in row]
        else:
            nonsingleton_sample_count = 0 
            nonsingleton_allele_count = 0 

        singleton_count = sum(c for b, c in fixed_counts.items() if c == 1)

    else:
        nonsingleton_sample_count = 0
        nonsingleton_allele_count = 0
        singleton_count = 0
    
    if singleton_count == 0:
        site_code = nonsingleton_allele_count
    elif nonsingleton_allele_count == 0:
        site_code = 0
    else:
        site_code = 5+nonsingleton_allele_count

    total_missing = missing_count + ploidy_fail_count

    return_count = {
        "Site_Code":site_code,
        "Missing":total_missing,
        "Uncovered":missing_count,
        "Filtered":ploidy_fail_count,
        "Fixed":nonsingleton_sample_count+singleton_count,
        "Hets":het_count,
        "Singletons":singleton_count,
        "Nonsingletons":nonsingleton_sample_count,
        "Nonsingleton_Alleles":nonsingleton_allele_count
    }
       
    return return_count, row, {total_missing:site_code}

def compile_chunks(sorted_files, out_file):
    for f in sorted_files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Missing chunk file: {f}")

    lazy_frames = [pl.scan_parquet(f) for f in sorted_files]
    (
        pl.concat(lazy_frames, how="vertical")
        .collect(streaming=True)
        .write_parquet(out_file, compression="snappy")
    )

    for f in sorted_files:
        os.remove(f)



##### MAIN #####
args = parse_args()

join_id = str(args.join_id)

# Output directory
output_directory = os.path.abspath(args.output_directory)
if not os.path.exists(output_directory):
    sys.exit(f"{output_directory} (--out_dir) does not exist")

temp_directory = os.path.join(output_directory,f"Temp_{join_id}")
if not os.path.exists(temp_directory):
    sys.exit(f"{temp_directory} does not exist")

# Scaffold file
scaffold_file = os.path.abspath(args.scaffold_file)
if not os.path.exists(scaffold_file):
    sys.exit(f"{scaffold_file} (--scaffold) does not exist")

# Base file
base_file = os.path.abspath(args.base_file)
if not os.path.exists(base_file):
    sys.exit(f"{base_file} (--base) does not exist")

output_site_file = os.path.join(output_directory,f"{join_id}_Sites.parquet")
output_code_file = os.path.join(output_directory,f"{join_id}_Codes.parquet")

if os.path.exists(output_site_file):
    sys.exit(f"{output_site_file} already exists...")
if os.path.exists(output_code_file):
    sys.exit(f"{output_code_file} already exists...")

# Create chunks
row_count = pq.ParquetFile(scaffold_file).metadata.num_rows

if args.mem_mode:
    max_chunks = os.cpu_count() * 40
else:
    max_chunks = os.cpu_count() * 4

n_chunks = min(max_chunks, row_count) if row_count > 0 else 1
chunk_size = (row_count + n_chunks - 1) // n_chunks if row_count > 0 else 0

jobs = []
for i in range(n_chunks):
    start = i * chunk_size
    stop = min(start + chunk_size, row_count)
    if start >= stop:
        continue

    jobs.append((i, start, stop, base_file, temp_directory))

code_list = []
site_list = []
dict_list = []

with ProcessPoolExecutor(max_workers=cpu) as ex:
    futures = {ex.submit(process_chunk, *job): job[0] for job in jobs}

    for fut in as_completed(futures):
        i = futures[fut]
        try:
            code_file, site_file, missing_dict = fut.result()
            code_list.append(code_file)
            site_list.append(site_file)
            dict_list.append(missing_dict)
        except Exception as e:
            print(f"[{i}] ❌ Error in worker: {e}")


# Save codes and sites
sorted_codes = natsorted(code_list)
sorted_sites = natsorted(site_list)

compile_chunks(sorted_codes,output_code_file)
compile_chunks(sorted_sites,output_site_file)

# Save missing summary

code_labels = {
    0: "Other",
    1: "Pure_Fixed",
    2: "Pure_Biallelic",
    3: "Pure_Triallelic",
    4: "Pure_Quadallelic",
    5: "Pure_Pentallelic",
    6: "Fixed_wSingleton",
    7: "Biallelic_wSingleton",
    8: "Triallelic_wSingleton",
    9: "Quadallelic_wSingleton",
}

missing_summary = defaultdict(Counter)
for chunk in dict_list:
    for m, codes in chunk.items():
        missing_summary[int(m)].update({int(k): int(v) for k, v in codes.items()})

df = (
    pl.DataFrame([
        {"Missing": m, **{code_labels.get(c, f"Code_{c}"): v for c, v in counts.items()}}
        for m, counts in missing_summary.items()
    ])
    .fill_null(0)
)

nonzero_cols = [
    "Missing",
    *[c for c in df.columns if c != "Missing" and df[c].sum() > 0]
]

missing_file = os.path.join(output_directory, f"{join_id}_Missing.tsv")
df.select(nonzero_cols).sort("Missing").write_csv(missing_file, separator="\t")

# region 04: Save Summary
site_count_file = os.path.join(output_directory,f"{join_id}_Site_Counts.tsv")
called_base_file = os.path.join(output_directory,f"{join_id}_Called_Bases.txt")
sample_ids = pl.scan_parquet(base_file).collect_schema().names()

joined_info = {
    "Join_ID":join_id,
    "Joined_Directory":output_directory,
    "Called_Base_File":called_base_file,
    "Sample_IDs":",".join(sample_ids),
    "Scaffold_File":scaffold_file,
    "Base_File":base_file,
    "Code_File":output_code_file,
    "Site_File":output_site_file,
    "Site_Count_File":site_count_file,
    "Missing_File":missing_file
    }

output_json_file = os.path.join(output_directory,f"{join_id}.json")

with open(output_json_file, "w", encoding="utf-8") as f:
    json.dump(joined_info, f, indent=4)