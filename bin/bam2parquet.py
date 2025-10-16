import polars as pl
import re
import os
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
from multiprocessing import Lock
import pyarrow.parquet as pq
import pyarrow as pa
from typing import List
import pandas as pd
import argparse
import json
import gzip
from datetime import datetime
import gc
from natsort import natsorted
from Bio import SeqIO
import subprocess
from io import StringIO
import pysam
import tempfile
import math
import csv
import io

lock = Lock()

def parse_args():
    parser = argparse.ArgumentParser(description="Process BAM/pileup file into SNPRS parquet")
    
    parser.add_argument("--bam", dest="bam_file", type=str, required=True,help="Path to input BAM file")
    parser.add_argument("--fasta", dest="fasta_file", type=str, required=True,help="Path to reference FASTA file")
    parser.add_argument("--parquet", dest="parquet_file", type=str, default=None,help="Path to output Parquet file [Default: bam_file/pileup file with .parquet extension]")
    
    parser.add_argument("--mapq", dest="mapq", type=int, default=15,help="Mapping quality argument for mpileup -q [Default: 15]")
    parser.add_argument("--baseq", dest="baseq", type=int, default=15,help="Base quality argument for mpileup -Q [Default: 15]")
    parser.add_argument("--adj_coef", dest="adj_coef", type=int, default=50,help="Adjustment coefficient for mpileup -C [Default: 50]")
    
    parser.add_argument("--cpus", dest="user_cpu", type=int, default=None,help="Number of CPUs requested [Default: all available CPUs]")

    return parser.parse_args()

def write_contig_bed(contigs, contig_lengths, bed_path):
    with open(bed_path, "w") as bed:
        for contig in contigs:
            length = contig_lengths.get(contig)
            if length is not None:
                bed.write(f"{contig}\t0\t{length}\n")
            else:
                sys.exit(f"ERROR: contig '{contig}' not found in contig_lengths")
                continue
    return bed_path


def run_mpileup_chunk(bam_file, fasta_file, mapq, baseq, adj_coef, bed_file, temp_file):
    
    cmd = (
        f"samtools mpileup -q {mapq} -Q {baseq} -C {adj_coef} -f {fasta_file} --ff 3844 --no-output-ends --no-output-del -l {bed_file} {bam_file} "
        f"| awk 'NF >= 4'"
    )

    results = [] 
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
        for line in proc.stdout:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue

            row = {
                "Scaffold": parts[0],
                "Position": int(parts[1]),
                "Ref_Base": parts[2],
                "Coverage": int(parts[3]),
                "Sample_Bases": parts[4],
            }

            for record in process_pileup_row(row):
                results.append([
                    record["contig_id"],
                    record["contig_position"],
                    record["depth"],
                    record["base"],
                    record["frequency"],
                ])

        stderr = proc.stderr.read()
        returncode = proc.wait()

    if returncode != 0:
        raise RuntimeError(f"mpileup failed for {bed_file}:\n{stderr}")

    os.remove(bed_file)

    if results:
        with lock, open(temp_file, "a", newline="") as out_f:
            writer = csv.writer(out_f, delimiter="\t")
            writer.writerows(results)

def process_pileup_row(row):
    contig = row["Scaffold"]
    pos = row["Position"]
    ref_base = row["Ref_Base"]
    sample_bases = row["Sample_Bases"]
    depth = row["Coverage"]

    freqs, total_depth = get_base_freqs(sample_bases, ref_base)

    for base, freq in freqs.items():
        yield {
            "contig_id": contig,
            "contig_position": pos,
            "depth": total_depth,
            "base": base,
            "frequency": freq,
        }

def get_base_freqs(sample_bases,ref_base):

    s = re.sub(r"-\d+", "", sample_bases)
    ins_counter = Counter()
    ins_pattern = r"\+(\d+)([ACGTNacgtn]+)"

    for match in re.finditer(ins_pattern, s):
        length = int(match.group(1))
        seq = match.group(2)[:length].upper()
        indel_str = f"+{length}{seq}"
        count = s.count(indel_str)
        ins_counter[indel_str] += count
 
    s = re.sub(ins_pattern, "", s).replace("*", "-").replace(".", ref_base).replace(",", ref_base).upper()
    depth = len(s)

    freqs = {b: c / depth for b, c in Counter(s).items()}
    indel_freqs = {indel: c / depth for indel, c in ins_counter.items()}
    freqs.update(indel_freqs)

    return freqs, depth

def compute_depth_stats(temp_file):

    depth_df = (
        pl.scan_csv(temp_file,separator="\t",has_header=True)
        .select(["contig_id", "contig_position", "depth"])
        .unique(subset=["contig_id", "contig_position"])
        .collect()
    )

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

    return stats

def compute_freq_stats(temp_file):

    freq_df = (
        pl.scan_csv(temp_file,separator="\t",has_header=True)
        .select(['frequency'])
        .collect()
    )
    
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
    
    return stats

def get_arrow(temp_file,contig_map,metadata):
    
    return (
        pl.scan_csv(temp_file,separator="\t",has_header=True)
        .with_columns([
            pl.col("contig_id").replace_strict(contig_map, default=None).alias("contig_index")
        ])
        .select(["contig_index", "contig_position", "base", "depth", "frequency"])
        .with_columns([
            pl.col("contig_index").cast(pl.Int64),
            pl.col("contig_position").cast(pl.Int64),
            pl.col("depth").cast(pl.Int64),
            pl.col("frequency").cast(pl.Float64),
            pl.col("base").cast(pl.Utf8)
        ])
        .sort(by=["contig_index", "contig_position", "frequency"],descending=[False, False, True])
        .collect()
        .to_arrow()
        .replace_schema_metadata(metadata)
    )

##### Main #####

# region 00: Parse args

args = parse_args()

user_cpu = args.user_cpu if args.user_cpu else os.cpu_count()
mapq = args.mapq
baseq = args.baseq
adj_coef = args.adj_coef

fasta_file = os.path.abspath(args.fasta_file)
if not os.path.exists(fasta_file):
    sys.exit(f"{fasta_file} does not exist...")

bam_file = os.path.abspath(args.bam_file)
if not os.path.exists(bam_file):
    sys.exit(f"{bam_file} does not exist...")
data_dir = os.path.dirname(bam_file)
sample_name = os.path.splitext(os.path.basename(bam_file))[0]

with pysam.AlignmentFile(bam_file, "rb") as bamfile:
    paired = any(read.is_paired for i, read in enumerate(bamfile.fetch(until_eof=True)) if i < 1000)

if not args.parquet_file:
    parquet_file = os.path.join(data_dir,sample_name+"_Raw.parquet")
else:
    parquet_file = os.path.abspath(args.parquet_file)

if os.path.exists(parquet_file):
    sys.exit(f"{parquet_file} already exists...")    

# endregion

# region 01: Process reference FASTA

raw_records = [(rec.id, str(rec.seq)) for rec in SeqIO.parse(fasta_file, "fasta")]

if not raw_records:
    sys.exit("No contigs found.")

raw_records = natsorted(raw_records, key=lambda x: x[0])    

total_sites = sum(len(seq) for _, seq in raw_records)

contig_ids = [rec_id.strip().split()[0] for rec_id, _ in raw_records]
contig_map = {v: k for k, v in enumerate(contig_ids)}
reverse_map = {k: v for k, v in enumerate(contig_ids)}

# endregion

# region 02: Process BAM into pileup

# Chunk out contigs by length
contig_lengths = {rec_id.strip().split()[0] : len(seq) for rec_id, seq in raw_records}

contig_count = len(contig_ids)
n_chunks = min(contig_count, user_cpu)
chunk_size = (contig_count + n_chunks - 1) // n_chunks

len_sorted_records = sorted(raw_records, key=lambda x: len(x[1]), reverse=True)
len_sorted_contig_ids = [rec_id.strip().split()[0] for rec_id, _ in len_sorted_records]
contig_chunks = [[] for _ in range(n_chunks)]

for i, contig in enumerate(len_sorted_contig_ids):
    contig_chunks[i % n_chunks].append(contig)

# Create BED files for mpileup
temp_beds = []
for i, chunk in enumerate(contig_chunks):
    bed_path = os.path.join(data_dir,sample_name+f"_tmp_{i}.bed")
    write_contig_bed(chunk, contig_lengths, bed_path)
    temp_beds.append(bed_path)

# Process pileup in parallel
temp_file = os.path.join(data_dir,sample_name+"_Temp.tsv")
with open(temp_file, "w", newline="") as out_f:
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow(["contig_id", "contig_position", "depth", "base", "frequency"])

with ProcessPoolExecutor(max_workers=user_cpu) as executor:
    futures = [
        executor.submit(
            run_mpileup_chunk,
            bam_file,
            fasta_file,
            mapq,
            baseq,
            adj_coef,
            bed_file,
            temp_file
        )
        for bed_file in temp_beds
    ]
    for fut in futures:
        try:
            fut.result()
        except Exception as e:
            print(f"❌ Error in worker: {e}")

with open(temp_file, "r") as f:
    line_count = sum(1 for _ in f)

if line_count == 0:
    raise ValueError("Pileup file was empty.")

# endregion

# region 04: Summarize and save 

depth_stats = compute_depth_stats(temp_file)
percent_covered = f"{int(depth_stats['covered'])/int(total_sites):.2f}"
freq_stats = compute_freq_stats(temp_file)

metadata = {
    "sample_id": sample_name,
    "bam_file": bam_file,
    "sample_parquet": parquet_file,
    "reference_genome": fasta_file,
    "percent_covered": percent_covered,
    "paired_end":"TRUE" if paired else "FALSE",
    "qc_filtering_params":f"MAPQ: {mapq}; BASEQ: {baseq}",
    "depth_statistics": ", ".join([f"{k}: {v}" for k, v in depth_stats.items()]),
    "allele_frequencies": ", ".join([f"{k}: {v}" for k, v in freq_stats.items()])
}

final_metadata = {k.encode(): str(v).encode() for k, v in metadata.items()}

try:
    arrow_table = get_arrow(temp_file, contig_map,final_metadata)
    pq.write_table(arrow_table, parquet_file, compression="snappy")
except Exception as e:
    raise e
else:
    os.remove(temp_file)

# endregion