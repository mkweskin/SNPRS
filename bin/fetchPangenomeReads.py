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
parser.add_argument('--read_dir',dest="read_dir",required=True, type=str, help='Path to directory containing nested read files')
parser.add_argument('--link_dir',dest="link_dir",default=None, type=str, help='Path to directory to link renamed files [Default: No links]')
parser.add_argument('--group',dest = "output",required=True, type=str, help='Path to group output csv')
parser.add_argument('--ext',dest="read_filetype",default='fastq.gz', type=str, help='Read extension')
parser.add_argument('--forward',dest = "forward_suffix",default='_1.fastq.gz', type=str, help='Forward suffix')
parser.add_argument('--reverse',dest = "reverse_suffix",default = '_2.fastq.gz', type=str, help='Reverse suffix')
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

if not paired_df["Sample_ID"].is_unique:
    seen = {}
    new_ids = []
    for sid in paired_df["Sample_ID"]:
        if sid not in seen:
            seen[sid] = 0
            new_ids.append(sid)
        else:
            seen[sid] += 1
            new_ids.append(f"{sid}_{seen[sid]}")
    paired_df["Sample_ID"] = new_ids
    
# Save group file
output_file = os.path.abspath(args.output)
paired_df.to_csv(output_file, sep=",", index=False,header=False)

# Link reads to link directory
if args.link_dir:
    
    link_dir = os.path.abspath(args.link_dir)
    
    for _, row in paired_df.iterrows():
        
        sample_id = row["Sample_ID"]
        
        group0 = row["Group_0"] if pd.notna(row["Group_0"]) else None
        group1 = row["Group_1"] if pd.notna(row["Group_1"]) else None
        
        fwd = row["Forward"]
        rev = row["Reverse"]

        name_parts = [sample_id] + [g for g in [group0, group1] if (g and g != "NO_GROUP")]
        base_name = "_".join(name_parts)

        fwd_name = f"{base_name}{forward_suffix}"
        rev_name = f"{base_name}{reverse_suffix}"
        se_name = f"{base_name}{read_filetype}"

        # Build full destination paths
        fwd_link = os.path.join(link_dir, fwd_name)
        rev_link = os.path.join(link_dir, rev_name)
        se_link = os.path.join(link_dir, se_name)

        if pd.isna(rev):
            se_link = os.path.join(link_dir, se_name)
            src_dst_pairs = [(fwd, se_link)]
        else:
            fwd_link = os.path.join(link_dir, fwd_name)
            rev_link = os.path.join(link_dir, rev_name)
            src_dst_pairs = [(fwd, fwd_link), (rev, rev_link)]

        for src, dst in src_dst_pairs:
            try:
                if os.path.islink(dst) or os.path.exists(dst):
                    os.remove(dst)
                os.symlink(os.path.abspath(src), dst)
            except Exception as e:
                print(f"Error linking {src} → {dst}: {e}")


# Pass [Sample_ID,Group_0,Group_1,Forward,Reverse] data out to Nextflow
paired_df.to_csv(sys.stdout, sep=",", index=False, header=False)

