#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import numpy as np

def group_reads(base_dir,read_filetype,forward_suffix,reverse_suffix,max_depth=3):
    records = []

    for root, dirs, files in os.walk(base_dir):
        rel_path = os.path.relpath(root, base_dir)
        depth = 0 if rel_path == "." else rel_path.count(os.sep) + 1

        if depth > max_depth:
            continue

        for f in files:
            if not f.endswith(read_filetype):
                continue

            rel_parts = [] if rel_path == "." else rel_path.split(os.sep)

            if not rel_parts:
                g0 = np.nan
            else:
                g0 = rel_parts[0]

            # Group_1
            if len(rel_parts) < 2:
                g1 = np.nan
            else:                    
                g1 = rel_parts[1]

            records.append({
                "Sample_ID": os.path.basename(f).replace(forward_suffix, "").replace(reverse_suffix, "").replace(read_filetype, ""),
                "Group_0": g0,
                "Group_1": g1,
                "Path": os.path.join(root, f)
            })

    return pd.DataFrame(records)

def pair_reads(group_df, read_filetype, forward_suffix, reverse_suffix):

    group_df = group_df.fillna("NO_GROUP")
    records = []

    def strip_suffix(path, suffixes):
        name = os.path.basename(path)
        for suf in suffixes:
            if name.endswith(suf):
                return name[: -len(suf)]
        return name

    for (g0, g1, sid), subdf in group_df.groupby(["Group_0", "Group_1", "Sample_ID"]):
        paths = list(subdf["Path"])

        fwd = [p for p in paths if p.endswith(forward_suffix)]
        rev = [p for p in paths if p.endswith(reverse_suffix)]
        others = [p for p in paths
                  if p.endswith(read_filetype)
                  and not p.endswith(forward_suffix)
                  and not p.endswith(reverse_suffix)]

        used_forwards = set()

        for r in rev:
            key = strip_suffix(r, [reverse_suffix])
            match = next((f for f in fwd if strip_suffix(f, [forward_suffix]) == key), None)

            if match:
                records.append({
                    "Sample_ID": sid,
                    "Group_0": g0,
                    "Group_1": g1,
                    "Forward": match,
                    "Reverse": r
                })
                used_forwards.add(match)
            else:
                records.append({
                    "Sample_ID": sid,
                    "Group_0": g0,
                    "Group_1": g1,
                    "Forward": r,
                    "Reverse": np.nan
                })

        for f in fwd:
            if f not in used_forwards:
                records.append({
                    "Sample_ID": sid,
                    "Group_0": g0,
                    "Group_1": g1,
                    "Forward": f,
                    "Reverse": np.nan
                })

        for o in others:
            records.append({
                "Sample_ID": sid,
                "Group_0": g0,
                "Group_1": g1,
                "Forward": o,
                "Reverse": np.nan
            })

    return pd.DataFrame(records)

# Parse args

parser = argparse.ArgumentParser(description='Fetch Reads')
parser.add_argument('-d','--dir',dest="read_dir", type=str, help='Path to directory containing nested read files')
parser.add_argument('-e','--extension',dest="read_filetype",default='fastq.gz', type=str, help='Read extension')
parser.add_argument('-f','--forward',dest = "forward_suffix",default='_1.fastq.gz', type=str, help='Forward suffix')
parser.add_argument('-r','--reverse',dest = "reverse_suffix",default = '_2.fastq.gz', type=str, help='Reverse suffix')
parser.add_argument('-o','--output',dest = "output",default = "NA", type=str, help='Path to group output csv')
args = parser.parse_args()

# Get read filetype information
read_filetype = args.read_filetype
if not read_filetype.startswith("."):
    read_filetype = "." + read_filetype

forward_suffix = args.forward_suffix
reverse_suffix = args.reverse_suffix

read_dir = os.path.abspath(args.read_dir)
if not os.path.exists(read_dir):
    sys.exit(f"Read directory {read_dir} does not exist...")

group_df = group_reads(read_dir, read_filetype,forward_suffix,reverse_suffix)
paired_df = pair_reads(group_df,read_filetype,forward_suffix,reverse_suffix)

if not args.output == "NA":
    paired_df.to_csv(args.output, sep=",", index=False,header=False)

paired_df.to_csv(sys.stdout, sep=",", index=False, header=False)