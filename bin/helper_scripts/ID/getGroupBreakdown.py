#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
import os
import numpy as np
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Estimate group mixture proportions from SNP data")
    parser.add_argument("--in", dest="snp_file", type=str, required=True, help="Path to SNP TSV file")
    parser.add_argument("--out", dest="out_file", type=str, required=True, help="Path to output TSV file")
    parser.add_argument("--error", dest="error", type=float, default=0.01, help="Small pseudo-count to avoid zeros")
    return parser.parse_args()


def load_all_samples_tsv(filename):

    data = defaultdict(dict)
    with open(filename, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row["Sample_ID"]
            gr = row["Reference_Species"]
            data[sample][gr] = {
                "snps": int(row["SNP_Count"]),
                "match": int(row["Match"]),
                "non_match": int(row["Non_Match"])
            }
    return dict(data)


def prepare_sample_data(group_data):

    groups = []
    match = []
    non_match = []

    for gr, counts in group_data.items():
        groups.append(gr)
        match.append(counts['match'])
        non_match.append(counts['non_match'])

    return groups, np.array(match, dtype=float), np.array(non_match, dtype=float)

def estimate_all_samples(all_samples, error=0.01, compute_ci=True):

    results = []

    for sample_id, group_data in all_samples.items():
        groups, match, non_match = prepare_sample_data(group_data)
        total = match + non_match

        p_raw = (match + error) / (total + 2 * error)
        p_hat = p_raw / np.sum(p_raw)

        if compute_ci:
            se_raw = np.sqrt(np.where(total > 0, p_raw * (1 - p_raw) / total, np.nan))
            se = se_raw / np.sum(p_raw)
            z = 1.96
            lower = np.maximum(0, p_hat - z * se)
            upper = np.minimum(1, p_hat + z * se)
        else:
            lower = upper = [None] * len(groups)

        for gr, p, lo, hi in zip(groups, p_hat, lower, upper):
            results.append({
                "Sample_ID": sample_id,
                "SNP_Group": gr,
                "Proportion": p,
                "CI_lower": lo,
                "CI_upper": hi
            })

    return pd.DataFrame(results)

args = parse_args()
snp_file = os.path.abspath(args.snp_file)
out_file = os.path.abspath(args.out_file)

all_samples = load_all_samples_tsv(snp_file)
df_results = estimate_all_samples(all_samples, error=args.error)

# Save to TSV
df_results.to_csv(out_file, sep="\t", index=False)


