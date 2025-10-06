#!/usr/bin/env python3

import os
import sys
import argparse
from pathlib import Path
import pysam

# Parse args

parser = argparse.ArgumentParser(description='Fetch Reads')
parser.add_argument('-b','--bam_data',dest="bam_data", type=str, help='Path to directory containing coordinate-sorted BAM files or a file with a list of BAM files')
args = parser.parse_args()

bam_data = os.path.abspath(args.bam_data)

bam_files = []
if os.path.isdir(bam_data):
    bam_files = [str(f.resolve()) for f in Path(bam_data).glob("*.bam")]
elif os.path.isfile(bam_data):
    with open(bam_data) as f:
        bam_files = [os.path.abspath(line.strip()) for line in f if line.strip() and line.strip().endswith(".bam")]

if len(bam_files) == 0:
    sys.exit(f"No BAM files detected via {bam_data}")

bam_tuples = [(os.path.splitext(os.path.basename(bam_file))[0], bam_file) for bam_file in bam_files]

for sample, bam in bam_tuples:
    with pysam.AlignmentFile(bam, "rb") as bamfile:
        sort_order = bamfile.header.get("HD", {}).get("SO", "unknown")
    if sort_order == "coordinate":        
        print(f"{sample},{bam}")
    else:
        print(f"Warning: {bam} is not coordinate-sorted and will not be processed", file=sys.stderr)