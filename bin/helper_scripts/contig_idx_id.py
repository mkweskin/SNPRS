import os
import sys
import argparse
from natsort import natsorted
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Process FASTA file and map contig IDs <-> indices")
    parser.add_argument("--fasta", dest="fasta_file", type=str, required=True, help="Path to reference FASTA file")
    parser.add_argument("--idx", dest="index", type=int, default=None, help="Contig index (integer)")
    parser.add_argument("--contig", dest="contig", type=str, default=None, help="Contig ID (string)")
    return parser.parse_args()


args = parse_args()

fasta_file = os.path.abspath(args.fasta_file)
if not os.path.exists(fasta_file):
    sys.exit(f"{fasta_file} does not exist...")

raw_records = [(rec.id.strip().split()[0], str(rec.seq)) for rec in SeqIO.parse(fasta_file, "fasta")]

if not raw_records:
    sys.exit("No contigs found.")

raw_records = natsorted(raw_records, key=lambda x: x[0])
contig_ids = [rec_id for rec_id, _ in raw_records]

contig_map = {v: k for k, v in enumerate(contig_ids)}
reverse_map = {k: v for k, v in enumerate(contig_ids)}

if args.index is not None:
    contig_id = reverse_map.get(args.index)
    if contig_id:
        print(f"Index {args.index} → Contig ID: {contig_id}")
    else:
        print(f"Could not find contig for index {args.index}")

if args.contig is not None:
    contig_idx = contig_map.get(args.contig)
    if contig_idx is not None:
        print(f"Contig ID {args.contig} → Index: {contig_idx}")
    else:
        print(f"Could not find index for contig {args.contig}")
