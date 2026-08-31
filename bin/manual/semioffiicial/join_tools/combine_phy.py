#!/usr/bin/env python3

import argparse
import glob
import os
import tempfile
import numpy as np
import polars as pl

def parse_args():
    p = argparse.ArgumentParser(
        description="Combine distance chunks into a PHYLIP distance matrix."
    )
    p.add_argument(
        "-i", "--input_folder",
        required=True,
        help="Folder containing dist_*.npy files"
    )
    p.add_argument(
        "-o", "--out",
        required=True,
        help="Output PHYLIP file"
    )
    p.add_argument(
        "-c", "--chunk_file",
        required=True,
        help="Output PHYLIP file"
    )
    return p.parse_args()


def main():

    args = parse_args()

    chunk_df = pl.read_csv(args.chunk_file, separator="\t")

    sample_ids = chunk_df["sample_id"].sort().to_list()
    chunk_ids = set(chunk_df["chunk_id"].sort().to_list())

    input_directory = os.path.abspath(args.input_folder)
    output_file = os.path.abspath(args.out)
    files = [os.path.join(input_directory,f"Chunk_Dist_{id}.npy") for id in chunk_ids]
    
    if not files:
        raise RuntimeError("No dist_*.npy files found.")
    
    N = len(sample_ids)

    first = np.load(files[0], mmap_mode="r")

    if first.shape[1] != N:
        raise RuntimeError("Chunk width does not match number of samples.")

    tmp = tempfile.NamedTemporaryFile(delete=False)
    tmp.close()

    matrix = np.memmap(
        tmp.name,
        dtype=np.float32,
        mode="w+",
        shape=(N, N)
    )

    matrix[:] = 0

    row = 0

    for f in files:

        block = np.load(f, mmap_mode="r")

        rows = block.shape[0]

        matrix[row:row + rows] = block

        row += rows

    if row != N:
        raise RuntimeError(
            f"Loaded {row:,} rows but expected {N:,}"
        )

    for i in range(N):
        matrix[i + 1:, i] = matrix[i, i + 1:]

    matrix.flush()

    with open(output_file, "w", buffering=1024 * 1024) as out:

        out.write(f"{N}\n")

        for i in range(N):

            out.write(sample_ids[i][:40].ljust(40))

            row = matrix[i]

            out.write(
                " ".join(f"{x:.6f}" for x in row)
            )

            out.write("\n")

    del matrix
    os.remove(tmp.name)
    for file in files:
        os.remove(file)
    print(output_file,end="")

if __name__ == "__main__":
    main()