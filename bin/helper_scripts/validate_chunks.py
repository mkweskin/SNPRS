import os
import sys
import pandas as pd

chunk_tsv = os.path.abspath(sys.argv[1])
chunk_df = pd.read_csv(chunk_tsv, sep="\t")

all_dir_files = {}

for _, row in chunk_df.iterrows():
    chunk_id = row["Chunk_ID"]
    chunk_directory = row["Chunk_Directory"]

    if not os.path.isdir(chunk_directory):
        raise AssertionError(f"Directory does not exist: {chunk_directory}")

    files = sorted(entry.name for entry in os.scandir(chunk_directory))
    all_dir_files[chunk_id] = set(files)

file_sets = list(all_dir_files.values())
first = file_sets[0]

for i, fset in enumerate(file_sets[1:], start=1):
    if fset != first:
        raise AssertionError(
            f"Mismatch in directory files:\n"
            f"Expected: {sorted(first)}\n"
            f"Found in {chunk_df.iloc[i]['Chunk_Directory']}: {sorted(fset)}"
        )

for _, row in chunk_df.iterrows():
    print(row["Chunk_Directory"])