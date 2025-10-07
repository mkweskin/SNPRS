import polars as pl
import re
import os
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
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

def parse_args():
    parser = argparse.ArgumentParser(description="Process BAM/pileup file into SNPRS parquet")
    
    parser.add_argument("--bam", dest="bam_file", type=str, required=True,help="Path to input BAM file")
    parser.add_argument("--pileup", dest="pileup_file", type=str, default=None,help="Optional: Path to input pileup file processed in SNPRS fashion")
    parser.add_argument("--parquet", dest="parquet_file", type=str, default=None,help="Path to output Parquet file [Default: bam_file/pileup file with .parquet extension]")
    parser.add_argument("--cpus", dest="user_cpu", type=int, default=None,help="Number of CPUs requested [Default: all available CPUs]")

    parser.add_argument("--fasta", dest="fasta_file", type=str, required=True,help="Path to reference FASTA file")

    parser.add_argument("--mapq", dest="mapq", type=int, default=15,help="Mapping quality argument for mpileup -q [Default: 15]")
    parser.add_argument("--baseq", dest="baseq", type=int, default=15,help="Base quality argument for mpileup -Q [Default: 15]")
    parser.add_argument("--adj_coef", dest="adj_coef", type=int, default=50,help="Adjusted coefficient argument for mpileup -C [Default: 50]")
    
    return parser.parse_args()

def process_chunk(path, i, start, end):
    
    data_dir = os.path.dirname(path)
    sample_name = os.path.splitext(os.path.basename(path))[0]

    chunk_file = os.path.join(data_dir,sample_name+f"_Chunk_{i}.parquet")

    columns = ["Scaffold", "Position", "Ref_Base", "Coverage", "Sample_Bases", "Code"]

    chunk_df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=columns,
        skiprows=start,
        nrows=end-start,
        dtype={
            "Scaffold": str,
            "Position": int,
            "Ref_Base": str,
            "Coverage": int,
            "Sample_Bases": str,
            "Code": str
        }
    ).drop(columns=["Code"])
    
    chunk_df = chunk_df[chunk_df["Coverage"] >= 1]

    records = []

    for _, row in chunk_df.iterrows():
        sample_bases = row["Sample_Bases"]
        ref_base = row["Ref_Base"]
        scaffold = row["Scaffold"]
        position = row["Position"]

        s = re.sub(r"-\d+", "", sample_bases)
        ins_counter = Counter()
        ins_pattern = r"\+(\d+)"
        
        while True:
            match = re.search(ins_pattern, s)
            if not match:
                break
            length = int(match.group(1))
            start_i = match.start()
            end_i = match.end() + length
            seq = s[match.end():end_i].upper()
            indel_str = f"+{length}{seq}"
            ins_counter[indel_str] += 1
            s = s[:start_i] + s[end_i:]
        
        s = s.replace("*", "-").replace(".", ref_base).replace(",", ref_base).upper()
        depth = len(s)
        freqs = {b: c / depth for b, c in Counter(s).items()}
        indel_freqs = {indel: c / depth for indel, c in ins_counter.items()}
        freqs.update(indel_freqs)

        base_rec = {
            "contig_id": scaffold.strip().split()[0],
            "contig_position": position,
            "depth": depth
        }
        
        for base, freq in freqs.items():
            rec = base_rec.copy()
            rec.update({"base": base, "frequency": freq})
            records.append(rec)

    col_types = {
        "contig_id": str,
        "contig_position": int,
        "base": str,
        "depth": int,
        "frequency": float
    }

    if not records:
        return None
    else:
        (
            pl.LazyFrame(records)
            .select(["contig_id", "contig_position", "base", "depth", "frequency"])
            .cast({
                "contig_id": pl.Utf8,
                "contig_position": pl.Int64,
                "base": pl.Utf8,
                "depth": pl.Int64,
                "frequency": pl.Float64,
            })
            .sink_parquet(chunk_file, compression="snappy")
        )

        return chunk_file

def convert_bam(bam_file,pileup_path,fasta_file,user_cpu, mapq, baseq, adj_coef):
    
    convert_cmd = f"samtools view -@ {user_cpu} -h -F 3844 {bam_file} | samtools mpileup -q {mapq} -Q {baseq} --no-output-ends --no-output-del -f {fasta_file} -C {adj_coef} - | awk 'NF >= 4 && ($4 > 0 || $4 != "")'"

    with open(pileup_path, "w") as pileup_file:
        subprocess.run(
            convert_cmd,
            shell=True,
            check=True,
            stdout=pileup_file,
            stderr=subprocess.DEVNULL,
            text=True
        )

def process_pileup(pileup_file,user_cpu,contig_map):
        
    with open(pileup_file, "r") as f:
        line_count = sum(1 for _ in f)

    if line_count == 0:
        os.remove(pileup_file)
        raise ValueError("Pileup file was empty.")

    n_chunks = min(line_count, user_cpu)
    chunk_size = (line_count + n_chunks - 1) // n_chunks
    
    jobs = []
    for i in range(n_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, line_count)
        jobs.append((pileup_file, i, start, end))
    
    # Process chunks
    chunk_files = []
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_chunk, path, i, start, end) for path, i, start, end in jobs]
        for future in as_completed(futures):
            chunk_files.append(future.result())

    return chunk_files
                        
def depth_stats(chunk_file):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/depth_stats.py")
    
    cmd = ["python", worker_script,chunk_file]
    result = subprocess.run(cmd,
        stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True
        )
    
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to process depth data.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
    try:
        output = json.loads(result.stdout)
    except json.JSONDecodeError:
        raise

    return output

def freq_stats(chunk_file):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    worker_script = os.path.join(script_dir, "helper_scripts/freq_stats.py")
    
    cmd = ["python", worker_script,chunk_file]
    result = subprocess.run(cmd,
        stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True
        )
    
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to process depth data.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return Code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
    try:
        output = json.loads(result.stdout)
    except json.JSONDecodeError:
        raise

    return output

def get_arrow(chunk_file,metadata):
    
    with open(chunk_file, "r") as f:
        files = [line.strip() for line in f if line.strip()]

    if not files:
        raise ValueError("No chunk files listed in " + chunk_file)

    lazy_frames = [pl.scan_parquet(f) for f in files]
    return (
        pl.concat(lazy_frames)
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
        .sort(['contig_index','contig_position'])
        .collect()
        .to_arrow()
        .replace_schema_metadata(metadata)
    )

#### MAIN ####
args = parse_args()

user_cpu = args.user_cpu if args.user_cpu else os.cpu_count()

mapq = args.mapq
baseq = args.baseq
adj_coef = args.adj_coef

# Process reference fasta
fasta_file = os.path.abspath(args.fasta_file)
if not os.path.exists(fasta_file):
    sys.exit(f"{fasta_file} does not exist...")
raw_records = [(rec.id, str(rec.seq)) for rec in SeqIO.parse(fasta_file, "fasta")]
if not raw_records:
    sys.exit("No contigs found.")
raw_records = natsorted(raw_records, key=lambda x: x[0])

total_sites = sum(len(seq) for _, seq in raw_records)
contig_ids = [rec_id.strip().split()[0] for rec_id, _ in raw_records]
contig_map = {v: k for k, v in enumerate(contig_ids)}

# Process BAM file
bam_file = os.path.abspath(args.bam_file)
if not os.path.exists(bam_file):
    sys.exit(f"{bam_file} does not exist...")
with pysam.AlignmentFile(bam_file, "rb") as bamfile:
    paired = any(read.is_paired for i, read in enumerate(bamfile.fetch(until_eof=True)) if i < 1000)

# Check for pileup file
if args.pileup_file:
    pileup_file = os.path.abspath(args.pileup_file)
    if not os.path.exists(pileup_file):
        sys.exit(f"{pileup_file} does not exist...")  
    data_dir = os.path.dirname(pileup_file)
    sample_name = os.path.splitext(os.path.basename(pileup_file))[0]
else:
    data_dir = os.path.dirname(bam_file)
    sample_name = os.path.splitext(os.path.basename(bam_file))[0]
    pileup_path = os.path.join(data_dir,sample_name+".pileup")
    pileup_file = convertBAM(bam_file,pileup_path,fasta_file,user_cpu,mapq,baseq,adj_coef)

# Set ouptut file
if not args.parquet_file:
    parquet_file = os.path.join(data_dir,sample_name+".parquet")
else:
    parquet_file = os.path.abspath(args.parquet_file)

if os.path.exists(parquet_file):
    sys.exit(f"{parquet_file} already exists...")    

# Convert pileup to chunks
chunk_files = process_pileup(pileup_file,user_cpu,contig_map)
chunk_file = os.path.join(data_dir,sample_name+"_Chunks.txt")

with open(chunk_file, mode="w+") as file:
    file.write("\n".join(chunk_files) + "\n")

depth_stats = depth_stats(chunk_file)
percent_covered = f"{int(depth_stats['covered'])/int(total_sites):.2f}"
freq_stats = freq_stats(chunk_file)

metadata = {
    "sample_id": sample_name,
    "bam_file": bam_file,
    "reference_genome": fasta_file,
    "percent_covered": percent_covered,
    "paired_end":"TRUE" if paired else "FALSE",
    "qc_filtering_params":f"MAPQ: {mapq}; BASEQ: {baseq}; ADJ_COEF: {adj_coef}",
    "depth_statistics": ", ".join([f"{k}: {v}" for k, v in depth_stats.items()]),
    "allele_frequencies": ", ".join([f"{k}: {v}" for k, v in freq_stats.items()])
}

final_metadata = {k.encode(): str(v).encode() for k, v in metadata.items()}

try:
    arrow_table = get_arrow(chunk_file, final_metadata)
    pq.write_table(arrow_table, parquet_file, compression="snappy")
except Exception as e:
    raise e
else:
    with open(chunk_file, "r") as f:
        files = [line.strip() for line in f if line.strip()]
    for f in files:
        try:
            os.remove(f)
        except FileNotFoundError:
            pass
    os.remove(chunk_file)




