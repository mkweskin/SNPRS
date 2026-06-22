#!/usr/bin/env python3

import os
import glob
import numpy as np
import argparse
import json
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-o", "--output_dir", default=".")
    parser.add_argument("-n", "--name", default="all_samples")
    parser.add_argument("-c", "--chunk", type=int, default=0)
    return parser.parse_args()


def main():

    args = parse_args()

    files = sorted(glob.glob(os.path.join(args.input_dir, "*.bin")))
    if not files:
        raise RuntimeError("No .bin files found")

    N = len(files)
    L = np.fromfile(files[0], dtype=np.int8).shape[0]

    chunk_ids = [i // args.chunk for i in range(N)] if args.chunk > 0 else [0]*N

    names = [os.path.splitext(os.path.basename(f))[0] for f in files]

    out_matrix  = os.path.join(os.path.abspath(args.output_dir), f"{args.name}_matrix.bin")
    out_names   = os.path.join(os.path.abspath(args.output_dir), f"{args.name}_names.tsv")
    out_json    = os.path.join(os.path.abspath(args.output_dir), f"{args.name}.json")

    M = np.memmap(out_matrix, dtype=np.int8, mode="w+", shape=(N, L))

    for i, f in enumerate(tqdm(files, desc="Loading samples", unit="sample")):
        M[i, :] = np.fromfile(f, dtype=np.int8)

    M.flush()
    del M

    with open(out_names, "w") as out:
        out.write("row_id\tchunk_id\tsample_id\n")
        for i, name in enumerate(names):
            out.write(f"{i}\t{chunk_ids[i]}\t{name}\n")

    with open(out_json, "w") as f:
        json.dump({
            "Sample_Count": N,
            "Alignment_Length": L,
            "Matrix_File": out_matrix,
            "Names_File": out_names
        }, f, indent=2)


if __name__ == "__main__":
    main()