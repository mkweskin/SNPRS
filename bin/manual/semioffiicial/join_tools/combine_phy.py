#!/usr/bin/env python3

import argparse
import glob
import os

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--indir",required=True,help="Directory containing dist_*.txt files")
    p.add_argument("--out",required=True,help="Output PHYLIP file")
    return p.parse_args()


def chunk_number(path):
    base = os.path.basename(path)
    return int(base.split("_")[1].split(".")[0])


def main():

    args = parse_args()

    files = sorted(
        glob.glob(os.path.join(args.indir, "dist_*")),
        key=chunk_number
    )

    if not files:
        raise RuntimeError("No chunk files found")

    print(f"[INFO] Found {len(files)} chunks")

    n_taxa = 0

    for f in files:
        with open(f) as fh:
            for line in fh:
                if line.strip():
                    n_taxa += 1

    print(f"[INFO] Total taxa: {n_taxa:,}")

    with open(args.out, "w", buffering=1024 * 1024) as out:

        out.write(f"{n_taxa}\n")

        written = 0

        for f in files:

            print(f"[INFO] Reading {f}")

            with open(f) as fh:

                for line in fh:

                    line = line.rstrip()

                    if not line:
                        continue

                    out.write(line)
                    out.write("\n")

                    written += 1

                    if written % 1000 == 0:
                        print(
                            f"[INFO] Wrote {written:,}/{n_taxa:,}",
                            flush=True
                        )

    print("[DONE]")


if __name__ == "__main__":
    main()
