import sys
import polars as pl
import os
import json
from collections import defaultdict, Counter
from natsort import natsorted
from concurrent.futures import ProcessPoolExecutor
import math
import argparse
import pyarrow.parquet as pq


def parse_args():
    parser = argparse.ArgumentParser(description="Create scaffold parquet from called base files")
    
    parser.add_argument("--scaffold", dest="scaffold_file", type=str, required=True,help="Path to scaffold parquet")
    parser.add_argument("--bases", dest="base_file", type=str, required=True,help="Path to base parquet")
    parser.add_argument("--out_dir", dest="output_directory", type=str, default=None,help="Path to output Parquet files [Default: cwd]")
    parser.add_argument("--join_id", dest="join_id", type=str, required=True,help="Prefix for output files")
    parser.add_argument("--mem_mode", dest="mem_mode", action = "store_true", help="Split chunks up 10X smaller than default")

    return parser.parse_args()

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

# Sample codes

# 0: Uncovered
# 1: Filtered
# 2: Het_Base
# 3: Het_Gap
# 4: Singleton_Base
# 5: Singleton_Gap
# 6: Nonsingleton_Base
# 7: Nonsingleton_Gap

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

def classify_row(row):
    fixed_bases = {'A', 'C', 'G', 'T', '-'}
    degenerate_codes = set("RYSWKMBDHV")
    het_bases = degenerate_codes
    gap_het_bases = set('actg') | set(c.lower() for c in degenerate_codes)

    fixed_counts = Counter(base for base in row if base in fixed_bases)
    singleton_bases = {b for b, c in fixed_counts.items() if c == 1}

    codes = []
    nonsingleton_alleles = set()

    for base in row:
        if base == "?":
            codes.append(0)
        elif base == "N":
            codes.append(1)
        elif base in het_bases:
            codes.append(2)
        elif base in gap_het_bases:
            codes.append(3)
        elif base in singleton_bases:
            codes.append(4 if base != "-" else 5)
        elif base in fixed_bases:
            nonsingleton_alleles.add(base)
            codes.append(6 if base != "-" else 7)
        else:
            codes.append(1)
    
    counts = Counter(codes)
    
    allele_count = len(nonsingleton_alleles)
    singleton_count = counts.get(4, 0) + counts.get(5, 0)
    missing_count = counts.get(0, 0) + counts.get(1, 0)

    if singleton_count == 0:
        site_code = allele_count
    elif allele_count == 0:
        site_code = 0
    else:
        site_code = 5+allele_count

    missing_dict = {missing_count:site_code}

    return_count = {
        "Site_Code":site_code,
        "Missing":missing_count,
        "Uncovered":counts.get(0, 0),
        "Filtered":counts.get(1, 0),
        "Hets":counts.get(2, 0) + counts.get(3, 0),
        "Het_Base":counts.get(2,0),
        "Het_Gap":counts.get(3,0),
        "Singletons":counts.get(4,0) + counts.get(5,0),
        "Singleton_Base":counts.get(4,0),
        "Singleton_Gap":counts.get(5,0),
        "Nonsingleton_Alleles":allele_count,
        "Nonsingletons":counts.get(6,0) + counts.get(7,0),
        "Nonsingleton_Base":counts.get(6,0),
        "Nonsingleton_Gap":counts.get(7,0)
    }
        
    return return_count, codes, missing_dict

def process_chunk(i, start, stop, base_file):
    
    data_dir = os.path.dirname(base_file)
    
    sample_ids = pl.scan_parquet(base_file).collect_schema().names()
    code_chunk_file = os.path.join(data_dir, f"Code_Chunk_{i}.parquet")
    site_chunk_file = os.path.join(data_dir, f"Site_Chunk_{i}.parquet")
    
    length = stop - start
    
    df_chunk = (
        pl.scan_parquet(base_file)
        .slice(start, length)
        .collect(streaming=True)
    )

    results = [classify_row(list(row)) for row in df_chunk.iter_rows()]
    counts_list, codes_list, missing_dicts = zip(*results)

    pl.DataFrame(counts_list, orient="row").write_parquet(site_chunk_file, compression="snappy")
    pl.DataFrame(codes_list, schema=sample_ids, orient="row").write_parquet(code_chunk_file, compression="snappy")

    # Aggregate missing dicts per chunk
    grouped_missing = defaultdict(Counter)
    for d in missing_dicts:
        for missing_count, code in d.items():
            grouped_missing[missing_count][code] += 1

    grouped_missing = {k: dict(v) for k, v in grouped_missing.items()}

    return code_chunk_file, site_chunk_file, grouped_missing

args = parse_args()

join_id = args.join_id

# Output directory
if args.output_directory is None:
    output_directory = os.getcwd()
else:
    output_directory = os.path.abspath(args.output_directory)
if not os.path.exists(output_directory):
    sys.exit(f"{output_directory} (--out_dir) does not exist")

# Scaffold file
scaffold_file = os.path.abspath(args.scaffold_file)
if not os.path.exists(scaffold_file):
    sys.exit(f"{scaffold_file} (--scaffold_file) does not exist")

# Base file
base_file = os.path.abspath(args.base_file)
if not os.path.exists(scaffold_file):
    sys.exit(f"{scaffold_file} (--base_file) does not exist")    

output_site_file = os.path.join(output_directory,f"{join_id}_Sites.parquet")
output_code_file = os.path.join(output_directory,f"{join_id}_Codes.parquet")

if os.path.exists(output_site_file):
    sys.exit(f"{output_site_file} already exists...")
if os.path.exists(output_code_file):
    sys.exit(f"{output_code_file} already exists...")

row_count = pq.ParquetFile(scaffold_file).metadata.num_rows

if params.mem_mode:
    n_chunks = min(row_count, (os.cpu_count()*40))
else:
    n_chunks = min(row_count, (os.cpu_count()*4))

chunk_size = (row_count + n_chunks - 1) // n_chunks

jobs = []
for i in range(n_chunks):
    start = i * chunk_size
    stop = min((i + 1) * chunk_size, row_count)
    jobs.append((i, start, stop, base_file))

code_list = []
site_list = []
dict_list = []

with ProcessPoolExecutor() as executor:
    futures = [executor.submit(process_chunk, i, start, stop, base_file) for i, start, stop, base_file in jobs]

    for fut in futures:
        try:
            code_file, site_file, missing_dict = fut.result()
            code_list.append(code_file)
            site_list.append(site_file)
            dict_list.append(missing_dict)
        except Exception as e:
            print(f"❌ Error in worker: {e}")

# Define code labels
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
    10: "Pentallelic_wSingleton"
}

# Save missing summary
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

# Save codes and sites
sorted_codes = natsorted(code_list)
sorted_sites = natsorted(site_list)

compile_chunks(sorted_codes,output_code_file)
compile_chunks(sorted_sites,output_site_file)