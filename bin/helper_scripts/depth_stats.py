#!/usr/bin/env python3
import argparse
import polars as pl
import sys
import json

def compute_depth_stats(chunk_file):

    with open(chunk_file, "r") as f:
        files = [line.strip() for line in f if line.strip()]

    if not files:
        raise ValueError("No chunk files listed in " + chunk_file)

    lazy_frames = [pl.scan_parquet(f) for f in files]
    
    depth_df = (
        pl.concat(lazy_frames)
        .select(["contig_id", "contig_position", "depth"])
        .unique(subset=["contig_id", "contig_position"])
    ).collect()

    covered_sites = depth_df.height

    q25 = round(depth_df["depth"].quantile(0.25, "nearest"), 2)
    q50 = round(depth_df["depth"].quantile(0.5, "nearest"), 2)
    q75 = round(depth_df["depth"].quantile(0.75, "nearest"), 2)
    mean_depth = round(depth_df["depth"].mean(), 2)
    min_depth = int(depth_df["depth"].min())
    max_depth = int(depth_df["depth"].max())

    stats = {
        "covered": str(covered_sites),
        "min": str(min_depth),
        "mean": str(mean_depth),
        "max": str(max_depth),
        "q25": str(q25),
        "q50": str(q50),
        "q75": str(q75)
    }

    print(json.dumps(stats))
    return stats

if __name__ == "__main__":
    chunk_file = sys.argv[1]
    compute_depth_stats(chunk_file)