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
    parser = argparse.ArgumentParser(description="Process pileup file into parquet")
    
    parser.add_argument("-b","--bam", dest="bam_file", type=str, required=True,help="Path to input BAM file")
    parser.add_argument("-f","--fasta", dest="fasta_file", type=str, required=True,help="Path to reference FASTA file")
    parser.add_argument("-p","--parquet", dest="parquet_file", type=str, default=None,help="Path to output Parquet file [Default: bam_file with .parquet extension]")

    parser.add_argument("--mapq", dest="mapq", type=int, default=15,help="Mapping quality argument for mpileup -q [Default: 15]")
    parser.add_argument("--baseq", dest="baseq", type=int, default=15,help="Base quality argument for mpileup -Q [Default: 15]")
    parser.add_argument("--adj_coef", dest="adj_coef", type=int, default=50,help="Adjusted coefficient argument for mpileup -C [Default: 50]")

    parser.add_argument("--dup", dest="duplicate", action = "store_true",help="Perform filtering of PCR duplicates")
    
    return parser.parse_args()

def process_chunk(path, start, end):
    
    pileup_schema = {
        "Scaffold": pl.Utf8,
        "Position": pl.Int64,
        "Ref_Base": pl.Utf8,
        "Coverage": pl.Int64,
        "Sample_Bases": pl.Utf8,
        "Code": pl.Utf8
    }

    chunk_df = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        schema_overrides=pileup_schema,
        skip_rows=start,
        n_rows=end - start
    ).drop("Code").filter(pl.col("Coverage") >= 1)

    records = []

    for row in chunk_df.iter_rows(named=True):
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

    if records:
        df = (pl.DataFrame(records)
            .select(["contig_id", "contig_position", "base", "depth", "frequency"])
            .cast({
                "contig_id": pl.Utf8,
                "contig_position": pl.Int64,
                "base": pl.Utf8,
                "depth": pl.Int64,
                "frequency": pl.Float64,
            })
        )        
        return df
    else:
        return pl.DataFrame(
            schema={
                "contig_id": pl.Utf8,
                "contig_position": pl.Int64,
                "base": pl.Utf8,
                "depth": pl.Int64,
                "frequency": pl.Float64,
            }
        )

def bam_to_pileup(bam_file,fasta_file,mapq,baseq,adj_coef,duplicate):
        
    cpus = os.cpu_count()

    bamfile = pysam.AlignmentFile(bam_file, "rb")
    bam_dir = os.path.dirname(os.path.abspath(bam_file))
    
    paired = any(read.is_paired for i, read in enumerate(bamfile.fetch(until_eof=True)) if i < 1000)
    bamfile.close()
    
    if duplicate:
        convert_cmd = f"samtools sort -n -@ {cpus} {bam_file} | samtools fixmate -@ {cpus} -m - - | samtools sort -@ {cpus} - | samtools markdup -@ {cpus} - - | samtools view -@ {cpus} -h -F 3844 - | samtools mpileup -q {mapq} -Q {baseq} --no-output-ends --no-output-del -f {fasta_file} -C {adj_coef} - | awk '$4 > 0'"
    else:
        convert_cmd = f"samtools view -@ {cpus} -h -F 3844 {bam_file} | samtools mpileup -q {mapq} -Q {baseq} --no-output-ends --no-output-del -f {fasta_file} -C {adj_coef} - | awk '$4 > 0'"

    with tempfile.NamedTemporaryFile(
        mode="w+", dir=bam_dir, suffix=".pileup", delete=False
    ) as tmp:
        subprocess.run(
            convert_cmd,
            shell=True,
            check=True,
            stdout=tmp,
            stderr=subprocess.DEVNULL,
            text=True
        )
        tmp_path = tmp.name

    with open(tmp_path, "r") as f:
        line_count = sum(1 for _ in f)

    if line_count == 0:
        os.remove(tmp_path)
        raise ValueError("Pileup file appears to be empty.")

    n_chunks = min(line_count, os.cpu_count())
    chunk_size = (line_count + n_chunks - 1) // n_chunks
    
    jobs = []
    for i in range(n_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, line_count)
        jobs.append((tmp_path, start, end))
    
    # Process chunks
    dfs = []
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_chunk, path, start, end) for path, start, end in jobs]
        for future in as_completed(futures):
            dfs.append(future.result())

    os.remove(tmp_path)
    return dfs, paired
                    
#### MAIN ####
args = parse_args()

# Check that files exist
fasta_file = os.path.abspath(args.fasta_file)
if not os.path.exists(fasta_file):
    sys.exit(f"{fasta_file} does not exist...")

bam_file = os.path.abspath(args.bam_file)
if not os.path.exists(bam_file):
    sys.exit(f"{bam_file} does not exist...")
bam_dir = os.path.dirname(bam_file)
sample_name = os.path.splitext(os.path.basename(bam_file))[0]

if not args.parquet_file:
    parquet_file = os.path.join(bam_dir,sample_name+".parquet")
else:
    parquet_file = os.path.abspath(args.parquet_file)

if os.path.exists(parquet_file):
    sys.exit(f"{parquet_file} already exists...")

# Process reference fasta
raw_records = [(rec.id, str(rec.seq)) for rec in SeqIO.parse(fasta_file, "fasta")]
if not raw_records:
    sys.exit("No contigs found.")
raw_records = natsorted(raw_records, key=lambda x: x[0])

total_sites = sum(len(seq) for _, seq in raw_records)
contig_ids = [rec_id.strip().split()[0] for rec_id, _ in raw_records]
contig_map = {v: k for k, v in enumerate(contig_ids)}

# Convert BAM to pileup
mapq = args.mapq
baseq = args.baseq
adj_coef = args.adj_coef

all_records,paired = bam_to_pileup(bam_file,fasta_file,mapq,baseq,adj_coef,args.duplicate)
out_df = (
    pl.concat(all_records)
    .with_columns([
        pl.col("contig_id").replace_strict(contig_map, default=None).alias("contig_index")
    ])
    .select(["contig_index","contig_position","base","depth","frequency"])
    .sort(['contig_index','contig_position'])
)

# Coverage info
depth_df = out_df.select(["contig_index", "contig_position","depth"]).unique()
overall_coverage = depth_df.height / total_sites

depth_col = depth_df["depth"]
depth_q25 = round(depth_col.quantile(0.25, "nearest"), 2)
depth_q50 = round(depth_col.quantile(0.5, "nearest"), 2)
depth_q75 = round(depth_col.quantile(0.75, "nearest"), 2)

major_allele_count = (out_df['frequency'] >= 0.85).sum()
intermediate_count = ((out_df['frequency'] >= 0.15) & (out_df['frequency'] < 0.85)).sum()
minor_allele_count = (out_df['frequency'] < 0.15).sum()

# Add metadata and write to Parquet
metadata = {
    "sample_id": sample_name,
    "bam_file": bam_file,
    "reference_genome": fasta_file,
    "percent_covered": str(overall_coverage),
    "paired_end":"TRUE" if paired else "FALSE",
    "duplicate_filtered": "TRUE" if args.duplicate else "FALSE",
    "qc_filtering_params":f"MAPQ: {mapq}; BASEQ: {baseq}; ADJ_COEF: {adj_coef}",
    "depth_quantiles": f"Q25: {depth_q25}; Q50: {depth_q50}; Q75: {depth_q75}",
    "allele_profiles": f"85: {major_allele_count}; 15: {intermediate_count}; Lt_15: {minor_allele_count}"
}

final_metadata = {k.encode(): v.encode() for k, v in metadata.items()}
out_df = out_df.to_arrow().replace_schema_metadata(final_metadata)
pq.write_table(out_df, parquet_file, compression="snappy")