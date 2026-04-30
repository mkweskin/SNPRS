import os
import sys
import argparse
from natsort import natsorted
from Bio import SeqIO
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Process FASTA file and map contig IDs <-> indices")
    
    # Read in inputs
    parser.add_argument("--fasta", dest="fasta_file", type=str, required=True, help="Path to reference FASTA file")
    parser.add_argument("--cpu_count", dest="cpu_count", type=int, required=True, help="CPU count (used to chunk data when converting BAM to parquet)")
    parser.add_argument("--bed_dir", dest="bed_dir", type=str, required=True, help="Output directory for BED files (must not exist)")

    return parser.parse_args()

def write_contig_bed(contigs, contig_lengths, bed_path, fai_path):

    with open(fai_path) as f:
        fai_order = [line.split('\t')[0] for line in f]

    order_index = {name: i for i, name in enumerate(fai_order)}

    try:
        sorted_contigs = sorted(contigs, key=lambda c: order_index[c])
    except KeyError as missing:
        sys.exit(f"ERROR: contig '{missing.args[0]}' not found in .fai file")

    with open(bed_path, "w") as bed:
        for contig in sorted_contigs:
            length = contig_lengths.get(contig)
            if length is None:
                sys.exit(f"ERROR: contig '{contig}' not found in contig_lengths")

            bed.write(f"{contig}\t0\t{length}\n")
            
args = parse_args()

fasta_file = os.path.abspath(args.fasta_file)
fai_file = fasta_file + ".fai"

if not os.path.exists(fasta_file):
    sys.exit(f"{fasta_file} does not exist...")
elif not os.path.exists(fai_file):
    sys.exit(f"{fai_file} does not exist...")

bed_dir = os.path.abspath(args.bed_dir)
if os.path.exists(bed_dir):
    sys.exit(f"{bed_dir} exists...")
try:
    os.mkdir(bed_dir)
except Exception:
    raise

fasta_name = os.path.splitext(os.path.basename(fasta_file))[0]
data_dir = os.path.dirname(fasta_file)

# Parse records
raw_records = [
    (rec.id.strip().split()[0], str(rec.seq))
    for rec in SeqIO.parse(fasta_file, "fasta")
]

if not raw_records:
    sys.exit("No contigs found.")

raw_records = natsorted(raw_records, key=lambda x: x[0])
contig_ids = [rec_id.strip().split()[0] for rec_id, _ in raw_records]
contig_map = {v: k for k, v in enumerate(contig_ids)}
reverse_map = {k: v for k, v in enumerate(contig_ids)}

contig_lengths = {
    rec_id.split()[0]: len(seq)
    for rec_id, seq in raw_records
}

len_sorted_records = sorted(raw_records, key=lambda x: len(x[1]), reverse=True)
len_sorted_contig_ids = [rec_id.split()[0] for rec_id, _ in len_sorted_records]

contig_count = len(contig_ids)
n_chunks = min(contig_count, args.cpu_count)
contig_chunks = [[] for _ in range(n_chunks)]

for i, contig in enumerate(len_sorted_contig_ids):
    contig_chunks[i % n_chunks].append(contig)

for i, chunk in enumerate(contig_chunks):
    bed_path = os.path.join(bed_dir,f"Chunk_{i}.bed")
    write_contig_bed(chunk, contig_lengths, bed_path, fai_file)