#!/usr/bin/env python3

import argparse
import json
import time
from multiprocessing import Pool
import numpy as np

# ============================================================
# CONFIG
# ============================================================

BLOCK = 4096

ALL = None
N = None
L = None


# ============================================================
# ARGUMENTS
# ============================================================

def parse_args():

    p = argparse.ArgumentParser()

    p.add_argument("--json", required=True)
    p.add_argument("--chunk-id", type=int, required=True)
    p.add_argument("--chunk-size", type=int, required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--processors",type=int,default=16)
    return p.parse_args()

def load_names(path):

    names = []

    with open(path) as f:

        next(f)

        for line in f:
            fields = line.rstrip().split("\t")
            names.append(fields[2])

    return (names)


def init_worker(path, shape):

    global ALL, N, L

    N, L = shape

    ALL = np.memmap(
        path,
        dtype=np.int8,
        mode="r",
        shape=(N, L)
    )

def compute_one(i):

    xi = ALL[i]

    dist = np.zeros(N, dtype=np.float32)
    cov = np.zeros(N, dtype=np.int32)

    for s in range(0, L, BLOCK):

        e = min(s + BLOCK, L)

        xib = xi[s:e]
        Xb = ALL[:, s:e]

        valid = (xib > 0) & (Xb > 0)

        cov += valid.sum(axis=1)

        diff = (Xb != xib) & valid

        dist += diff.sum(axis=1)

    dist /= np.maximum(cov, 1)

    dist[i] = 0.0

    return i, dist

def main():

    args = parse_args()

    with open(args.json) as f:
        meta = json.load(f)

    N = meta["Sample_Count"]
    L = meta["Alignment_Length"]

    names = load_names(meta["Names_File"])

    start = args.chunk_id * args.chunk_size
    end = min(start + args.chunk_size, N)

    total = end - start

    print(
        f"[INFO] Chunk {args.chunk_id}: "
        f"samples {start:,} - {end-1:,} "
        f"({total:,} isolates)",
        flush=True
    )

    print(
        f"[INFO] Matrix shape: "
        f"{N:,} x {L:,}",
        flush=True
    )

    print(
        f"[INFO] CPUs: {args.processors}",
        flush=True
    )

    t0 = time.time()

    completed = 0

    with open(args.out, "w", buffering=1024 * 1024) as out:

        with Pool(
            processes=args.processors,
            initializer=init_worker,
            initargs=(
                meta["Matrix_File"],
                (N, L)
            ),
            maxtasksperchild=500
        ) as pool:

            for i, dist in pool.imap(
                compute_one,
                range(start, end),
                chunksize=16
            ):

                out.write(
                    names[i].ljust(40)[:40]
                )

                out.write(
                    " ".join(
                        f"{x:.6f}"
                        for x in dist
                    )
                )

                out.write("\n")

                completed += 1

                if (
                    completed % 25 == 0
                    or completed == total
                ):

                    elapsed = time.time() - t0

                    rate = completed / max(elapsed, 1)

                    remaining = total - completed

                    eta_seconds = remaining / max(rate, 1e-9)

                    eta_hours = eta_seconds / 3600

                    pct = (
                        completed
                        / total
                        * 100
                    )

                    print(
                        f"[Chunk {args.chunk_id}] "
                        f"{completed:,}/{total:,} "
                        f"({pct:.1f}%) | "
                        f"{rate:.2f} isolates/sec | "
                        f"ETA {eta_hours:.2f} h",
                        flush=True
                    )

    elapsed = time.time() - t0

    print(
        f"[DONE] Chunk {args.chunk_id} "
        f"completed in "
        f"{elapsed/3600:.2f} h",
        flush=True
    )


if __name__ == "__main__":
    main()
