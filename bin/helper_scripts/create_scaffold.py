import sys
import polars as pl
import os
import json
from natsort import natsorted
import pyarrow.parquet as pq
import psutil

def print_memory_usage():
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 * 1024)  # in MiB
    print(f"🔍 Memory usage: {mem:.2f} MiB")
    
def fetch_base_parquets(file_path):
    if not os.path.isfile(file_path):
        print(f"Error: Input file '{file_path}' does not exist.")
        sys.exit(1)

    with open(file_path, "r") as f:
        paths = [line.strip() for line in f if line.strip()]

    for path in paths:
        if not os.path.exists(path):
            print(f"Error: Parquet file '{path}' not found.")
            sys.exit(1)
    if len(paths) < 2:
        print("Error: You must provide at least two called base parquet files.")
        sys.exit(1)

    return [os.path.abspath(path) for path in paths]

def write_scaffold_parquet(output_directory,temp_directory,sorted_samples,sorted_parquet_paths,analysis_name):
    
    scaffold_parquet = os.path.join(output_directory,f"{analysis_name}_Scaffold.parquet")
    chunk_json = os.path.join(temp_directory,"Full_Chunks.json")
    
    scaffold = ( pl.scan_parquet(sorted_parquet_paths[0]) .filter(pl.col("type").is_in([0, 1, 4])) .select(["contig_index", "contig_position"]) .unique() .collect() ) 
    
    for path in sorted_parquet_paths[1:]: 
        df = ( pl.scan_parquet(path) .filter(pl.col("type").is_in([0, 1, 4])) 
              .select(["contig_index", "contig_position"]) 
              .unique() 
              .collect() 
              ) 
        new_rows = df.join(scaffold, on=["contig_index", "contig_position"], how="anti") 
        if new_rows.height > 0: 
            scaffold = pl.concat([scaffold, new_rows]) 
                
    scaffold.sort(['contig_index','contig_position']).write_parquet(scaffold_parquet)

    site_count = scaffold.height
    n_chunks = min(site_count, os.cpu_count()*4)
    chunk_size = max(1, site_count // n_chunks)
    remainder = site_count % n_chunks

    chunk_coords = []
    start = 0
    for i in range(n_chunks):
        size = chunk_size + (1 if i < remainder else 0)
        end = start + size
        chunk_coords.append((start, end))
        start = end

    if len(chunk_coords) >= 2:
        first_size = chunk_coords[0][1] - chunk_coords[0][0]
        last_size = chunk_coords[-1][1] - chunk_coords[-1][0]
        if last_size < 0.25 * first_size:
            prev_start, _ = chunk_coords[-2]
            _, last_end = chunk_coords[-1]
            chunk_coords = chunk_coords[:-2] + [(prev_start, last_end)]

    chunk_indexes = list(enumerate(chunk_coords))

    with open(chunk_json, "w") as f:
        json.dump(chunk_indexes, f)
        
    output = {
    "scaffold_parquet":scaffold_parquet,
    "chunk_indexes": chunk_json,
    "sorted_samples":sorted_samples,
    "sorted_parquets":sorted_parquet_paths
    }

    print(json.dumps(output)) 

if __name__ == "__main__":
    parquet_path_file,output_directory,temp_directory,analysis_name = sys.argv[1], sys.argv[2],sys.argv[3],sys.argv[4]
    
    parquet_files = fetch_base_parquets(parquet_path_file)

    sample_id_paths = [(os.path.basename(path).replace("_Called.parquet", ""), path) for path in parquet_files]
    sorted_path_sample_pairs = natsorted(sample_id_paths, key=lambda x: x[0])
    sorted_samples, sorted_parquet_paths = zip(*sorted_path_sample_pairs)
    
    write_scaffold_parquet(output_directory,temp_directory,sorted_samples,sorted_parquet_paths,analysis_name)

