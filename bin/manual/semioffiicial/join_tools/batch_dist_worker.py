#!/usr/bin/env python3


import argparse
import json
import numpy as np
from multiprocessing import Pool
import os
import polars as pl

SAMPLE_BLOCK = 64
QUERY_BLOCK = 8

MATRIX = None
N = None
L = None


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--chunk_file", required=True)
    p.add_argument("--chunk_size", type=int, required=True)
    p.add_argument("--chunk_id", type=int, required=True)

    p.add_argument("--scaffold", required=True)
    p.add_argument("--matrix", required=True)


    p.add_argument("--out", required=True)
    p.add_argument("--out_dir", required=True)

    p.add_argument("--raw", action="store_true")
    p.add_argument("--processors", type=int, default=1)
    return p.parse_args()


def load_names(path):
    names = []
    with open(path) as f:
        next(f)
        for line in f:
            names.append(line.rstrip().split("\t")[2])
    return names

def init_worker(path, shape):
    global MATRIX, N, L
    N, L = shape
    MATRIX = np.memmap(
        path,
        dtype=np.int8,
        mode="r",
        shape=(N, L)
    )

def compute_chunk(start, end, raw):

    qsize = end - start

    distances = np.zeros(
        (qsize, N),
        dtype=np.float32
    )

    coverage = np.zeros(
        (qsize, N),
        dtype=np.int32
    )

    query = MATRIX[start:end]

    for s in range(start, N, SAMPLE_BLOCK):

        e = min(s + SAMPLE_BLOCK, N)

        block = MATRIX[s:e]

        valid = (
            (query[:, None, :] > 0)
            &
            (block[None, :, :] > 0)
        )

        cov = valid.sum(axis=2)

        diff = (
            (query[:, None, :] != block[None, :, :])
            &
            valid
        )

        snps = diff.sum(axis=2)

        for q in range(qsize):

            global_q = start + q

            offset_start = max(s, global_q + 1)

            if offset_start < e:

                local_start = offset_start - s

                distances[
                    q,
                    offset_start:e
                ] += snps[
                    q,
                    local_start:e-s
                ]

                coverage[
                    q,
                    offset_start:e
                ] += cov[
                    q,
                    local_start:e-s
                ]

    if not raw:
        distances /= np.maximum(coverage, 1)

    for i in range(qsize):
        distances[i, start+i] = 0

    return start, distances

def run_task(args):
    return compute_chunk(*args)

def main():

    args = parse_args()

    chunk_file = os.path.abspath(args.chunk_file)
    scaffold_file = os.path.abspath(args.scaffold)
    matrix_file = os.path.abspath(args.matrix)
    chunk_id = int(args.chunk_id)
    chunk_size = int(args.chunk_size)
    
    data_dir = os.path.abspath(os.path.dirname(matrix_file))
    chunk_df = pl.read_csv(chunk_file, separator="\t")

    N = chunk_df.height
    
    L = (
        pl.scan_parquet(scaffold_file)
        .select(pl.len())
        .collect()
        .item()
    )

    chunk_start = chunk_id * chunk_size
    chunk_end = min(chunk_start + chunk_size, N)

    tasks = []

    for s in range(chunk_start, chunk_end, QUERY_BLOCK):

        e = min(s + QUERY_BLOCK, chunk_end)

        tasks.append(
            (s, e, args.raw)
        )

    results = {}

    with Pool(
        processes=args.processors,
        initializer=init_worker,
        initargs=(
            matrix_file,
            (N, L)
        ),
        maxtasksperchild=10
    ) as pool:

        for row_start, dist in pool.imap_unordered(
            run_task,
            tasks
        ):

            results[row_start] = dist
            
    blocks = [
        results[row_start]
        for row_start in sorted(results)
    ]

    block = np.vstack(blocks)

    np.save(
        args.out,
        block
    )
    
    print(data_dir,end="")
        
if __name__ == "__main__":
    main()