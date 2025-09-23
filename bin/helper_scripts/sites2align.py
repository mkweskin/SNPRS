import polars as pl
import json
import sys

def make_alignment(base_file,row_file,output_fasta,chunk_size=100000):

    with open(row_file, "r", encoding="utf-8") as f:
        surviving_rows = set(json.load(f))

    bases = pl.scan_parquet(base_file)
    sample_ids = bases.collect_schema().names()
    bases_filtered = bases.with_row_index("row_nr").filter(pl.col("row_nr").is_in(list(surviving_rows))).collect().drop("row_nr")

    # Stream FASTA
    with open(output_fasta, "w", encoding="utf-8") as f:
        for batch in bases_filtered.iter_slices(n_rows=chunk_size):
            for sid in sample_ids:
                seq_chunk = batch.get_column(sid).to_list()
                if seq_chunk:
                    f.write(f">{sid}\n{''.join(seq_chunk)}\n")
                    
if __name__ == "__main__":

    base_file,row_file,output_fasta= sys.argv[1],sys.argv[2],sys.argv[3]
    make_alignment(base_file,row_file,output_fasta)