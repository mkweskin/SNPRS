
#!/usr/bin/env python3

import polars as pl
import os
import sys
import argparse
import json
from Bio import SeqIO,AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from natsort import natsorted

parser = argparse.ArgumentParser(description='Generate alignment from SNPRS data')
parser.add_argument('-j','--join_dir',dest="join_dir", type=str,required=True, help='Path to SNPRS joined directory')
parser.add_argument('-m','--missing',dest="missing", type=str,required=True, help='If an integer, the max number allowed missing. If < 1, the proportion with data required')
parser.add_argument('-a','--alignment',dest="alignment_file", required=True,type=str, help='Output path for alignment')
parser.add_argument('-g','--gaps',dest="use_gaps",action="store_true", help='Include gaps as - in alignment, and use "N" for everything else (missing or invalid)')

def interpret_missing(missing_str, sample_count):

    try:
        val = int(missing_str)
        return val
    except ValueError:
        val = float(missing_str)

        if 0 < val < 1:
            required = int(sample_count * val)
            allowed_missing = sample_count - required
            return allowed_missing
        else:
            raise ValueError(
                f"--missing value '{missing_str}' must be an integer or a proportion between 0 and 1."
            )
            
args = parser.parse_args()

join_dir = os.path.abspath(args.join_dir)
join_id = os.path.basename(join_dir)

scaffold_file = os.path.join(join_dir,join_id,"_Scaffold.parquet")
code_file =  os.path.join(join_dir,f"{join_id}_Codes.parquet")
site_file =  os.path.join(join_dir,f"{join_id}_Sites.parquet")

lazy_codes = pl.scan_parquet(code_file)
lazy_sites = pl.scan_parquet(site_file)

cols = lazy_codes.collect_schema().names()
sample_ids = natsorted([c for c in cols if c not in {"contig_index", "contig_position"}])
sample_count = len(sample_ids)
max_missing = interpret_missing(args.missing,sample_count)

# Create Alignment
tree_sites = (
    lazy_sites
    .filter(pl.col("Nonsingleton_Alleles") >= 2)
    .filter(pl.col("Filtered") == 0)
    .filter(pl.col("Missing") <= max_missing)
    .select(["contig_index", "contig_position", "Missing"])
)

base_convert_dict = { 0:'N',
                    1:'A',2:"C",3:"G",4:"T",16:"-",
                    33:'A', 34:'C', 35:'T', 36:'G'}

if args.use_gaps:
    base_set = {0,1,2,3,4,16,33,34,35,36}
else:
    base_set = {0,1,2,3,4,33,34,35,36}

allowed = pl.lit(list(base_set))

lazy_tree_base_codes = (
    tree_sites
    .join(lazy_codes, on=["contig_index","contig_position"], how="left")
    .filter(
        pl.fold(
            acc=True,
            exprs=[pl.col(s).is_in(allowed) for s in sample_ids],
            function=lambda acc, x: acc & x
        )
    )
).collect(engine="streaming")

base_merged_df = lazy_tree_base_codes.with_columns([
    pl.col(col)
      .replace_strict(base_convert_dict)   
      .alias(col)
    for col in sample_ids
])

base_seq_records = []
for sample in sample_ids:
    sequence = "".join(base_merged_df[sample].to_list()) 
    seq_record = SeqRecord(Seq(sequence), id=sample)
    base_seq_records.append(seq_record)

base_alignment = MultipleSeqAlignment(base_seq_records)

alignment_file = os.path.abspath(args.alignment_file)
AlignIO.write(base_alignment, alignment_file, "fasta")