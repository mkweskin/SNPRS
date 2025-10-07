#!/usr/bin/env python3
import argparse
import polars as pl
import sys
import json

def computer_freq_stats(chunk_file):

    with open(chunk_file, "r") as f:
        files = [line.strip() for line in f if line.strip()]

    if not files:
        raise ValueError("No chunk files listed in " + chunk_file)

    lazy_frames = [pl.scan_parquet(f) for f in files]
    
    freq_df = (
        pl.concat(lazy_frames)
        .select(["frequency"])
    ).collect()

    freqs = freq_df["frequency"]

    total_alleles = freqs.len()

    stats = {
        "Allele_Count": str(total_alleles),
        "BT_0_1": f"{(freqs < 0.01).sum()}",
        "BT_1_5": f"{((freqs >= 0.01) & (freqs < 0.05)).sum()}",
        "BT_5_10": f"{((freqs >= 0.05) & (freqs < 0.10)).sum()}",
        "BT_10_15": f"{((freqs >= 0.10) & (freqs < 0.15)).sum()}",
        "BT_15_85": f"{((freqs >= 0.15) & (freqs < 0.85)).sum()}",
        "BT_85_90": f"{((freqs >= 0.85) & (freqs < 0.90)).sum()}",
        "BT_90_95": f"{((freqs >= 0.90) & (freqs < 0.95)).sum()}",
        "BT_95_99": f"{((freqs >= 0.95) & (freqs < 0.99)).sum()}",
        "BT_99_100": f"{(freqs >= 0.99).sum()}",
    }
    
    print(json.dumps(stats))
    return stats

if __name__ == "__main__":
    chunk_file = sys.argv[1]
    computer_freq_stats(chunk_file)