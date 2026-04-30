import os
import sys
import argparse
from natsort import natsorted
from Bio import SeqIO
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Process FASTA file and map contig IDs <-> indices")
    
    # Read in inputs
    parser.add_argument("--fasta", dest="fasta_file", type=str, default=None, help="Path to reference FASTA file")
    parser.add_argument("--parquet", dest="parquet_file", type=str, default=None, help="Path to reference parquet file")
    parser.add_argument("--cpu_count", dest="cpu_count", type=int, default=1, help="CPU count (used to chunk data when converting BAM to parquet; Default: 1)")

    # Generate output
    parser.add_argument("--make_parquet", dest="make_parquet", action="store_true", help="Create reference parquet with contig ID, index, and length from FASTA (requires --fasta)")
    parser.add_argument("--make_tsv", dest="make_tsv", action="store_true", help="Create reference TSV with contig ID, index, and length (requires --fasta or --parquet)")

    # Convert ID/IDX
    parser.add_argument("--idx", dest="index", type=int, default=None, help="Contig index (integer) or a file with a list of contig indexes")
    parser.add_argument("--contig", dest="contig", type=str, default=None, help="Contig ID (string) or a file with a list of contig IDs")

    return parser.parse_args()

def read_input_values(request):
    if os.path.isfile(request):
        with open(request) as f:
            return [line.strip() for line in f if line.strip()]
    else:
        return [request]

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

# Check run modes
if (args.index or args.contig) and (args.make_parquet or args.make_tsv):
    sys.exit("Script can be run in either fetch mode (--contig/--idx) or generate mode (--make_parquet/--make_tsv), but not both.")

if args.make_parquet and not args.fasta_file:
    sys.exit("If running in generate mode, reference information must be provided via --fasta")

# Process FASTA file
if args.fasta_file:  
    
    fasta_file = os.path.abspath(args.fasta_file)
    fai_file = fasta_file + ".fai"
    
    if not os.path.exists(fasta_file):
        sys.exit(f"{fasta_file} does not exist...")
    elif not os.path.exists(fai_file):
        sys.exit(f"{fai_file} does not exist...")

    fasta_name = os.path.splitext(os.path.basename(fasta_file))[0]
    data_dir = os.path.dirname(fasta_file)

    # Parse records
    raw_records = [
        (rec.id.strip().split()[0], str(rec.seq))
        for rec in SeqIO.parse(fasta_file, "fasta")
    ]

    if not raw_records:
        sys.exit("No contigs found.")

    # Sort contigs naturally (e.g. contig1, contig2, ..., contig10)
    raw_records = natsorted(raw_records, key=lambda x: x[0])
    contig_ids = [rec_id.strip().split()[0] for rec_id, _ in raw_records]
    contig_map = {v: k for k, v in enumerate(contig_ids)}
    reverse_map = {k: v for k, v in enumerate(contig_ids)}
    
    # Create DataFrame matching Parquet format
    df = pd.DataFrame({
        "contig_idx": range(len(contig_ids)),
        "contig_id": contig_ids,
        "contig_length": [len(seq) for _, seq in raw_records],
    })

# Process input Parquet file
elif args.parquet_file:
    parquet_file = os.path.abspath(args.parquet_file)

    if not os.path.exists(parquet_file):
        sys.exit(f"{parquet_file} does not exist...")

    fasta_name = os.path.splitext(os.path.basename(parquet_file))[0]
    data_dir = os.path.dirname(parquet_file)

    df = pd.read_parquet(
        parquet_file,
        columns=["contig_idx", "contig_id", "contig_length"]
    )
    contig_map = dict(zip(df["contig_id"], df["contig_idx"]))
    reverse_map = dict(zip(df["contig_idx"], df["contig_id"]))

else:
    sys.exit("You must provide either --fasta_file or --parquet_file.")

# If making the parquet, also make the chunked BED files for mpileup
if args.make_parquet:
    
    out_parquet = os.path.join(data_dir,fasta_name+".parquet")
    bed_dir = os.path.join(data_dir,"BED_Chunks")

    if os.path.exists(out_parquet):
        sys.exit(f"Output file {out_parquet} exists...")
    elif os.path.exists(bed_dir):
        sys.exit(f"BED directory {bed_dir} exists...")
    else:
        df.to_parquet(out_parquet, index=False,compression = "snappy")
        os.mkdir(bed_dir)

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

if args.make_tsv:
    out_tsv = os.path.join(data_dir,fasta_name+".tsv")
    df.to_csv(out_tsv, sep="\t", index=False)

if args.index is not None:
    index_values = read_input_values(str(args.index))
    for idx in index_values:
        try:
            idx_int = int(idx)
        except ValueError:
            print(f"Skipping invalid index: {idx}")
            continue

        contig_id = reverse_map.get(idx_int)
        if contig_id:
            print(f"Index {idx_int} → Contig ID: {contig_id}")
        else:
            print(f"Could not find contig for index {idx_int}")

if args.contig is not None:
    contig_values = read_input_values(str(args.contig))
    for contig in contig_values:
        contig_idx = contig_map.get(contig)
        if contig_idx is not None:
            print(f"Contig ID {contig} → Index: {contig_idx}")
        else:
            print(f"Could not find index for contig {contig}")