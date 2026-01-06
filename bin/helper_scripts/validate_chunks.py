import os
import sys
import pandas as pd
import pyarrow.parquet as pq
from collections import defaultdict

chunk_tsv = os.path.abspath(sys.argv[1])
chunk_df = pd.read_csv(chunk_tsv, sep="\t")

all_dir_files = {}
parquet_row_counts = defaultdict(dict)

for _, row in chunk_df.iterrows():
    chunk_id = row["Chunk_ID"]
    chunk_directory = row["Chunk_Directory"]

    if not os.path.isdir(chunk_directory):
        raise AssertionError(f"Directory does not exist: {chunk_directory}")

    files = sorted(os.listdir(chunk_directory))
    all_dir_files[chunk_id] = set(files)

    for fname in files:
        if fname.endswith(".parquet"):
            fpath = os.path.join(chunk_directory, fname)
            try:
                n_rows = pq.ParquetFile(fpath).metadata.num_rows
            except Exception as e:
                raise AssertionError(f"Failed reading {fpath}: {e}")

            parquet_row_counts[chunk_id][fname] = n_rows

file_sets = list(all_dir_files.values())
first = file_sets[0]
for i, fset in enumerate(file_sets[1:], start=1):
    if fset != first:
        raise AssertionError(
            f"Mismatch in directory files:\n"
            f"Expected: {sorted(first)}\n"
            f"Found in {chunk_df.iloc[i]['Chunk_Directory']}: {sorted(fset)}"
        )

for chunk_id, counts in parquet_row_counts.items():
    uniq = set(counts.values())
    if len(uniq) > 1:
        raise AssertionError(
            f"Parquet row mismatch within chunk {chunk_id}:\n"
            + "\n".join([f"  {fname}: {n}" for fname, n in counts.items()])
        )
        
for _, row in chunk_df.iterrows():
    print(row["Chunk_Directory"])
