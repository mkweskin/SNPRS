import polars as pl
import sys
import textwrap
import os

def append_alignment(og_alignment_file, og_scaffold_file, called_base_files,new_alignment_file):

    with open(called_base_files) as f:
        called_files = [line.strip() for line in f if line.strip()]

    missing_files = [f for f in called_files if not os.path.exists(f)]
    if missing_files:
        raise FileNotFoundError(f"The following files are missing: {', '.join(missing_files)}")

    called_ids = [
        os.path.basename(f).replace("_Called.parquet", "")
        for f in called_files
    ]

    fasta_ids = []
    with open(og_alignment_file) as f:
        for line in f:
            if line.startswith(">"):
                fasta_ids.append(line[1:].strip().split()[0])

    overlap = set(called_ids) & set(fasta_ids)
    if overlap:
        raise ValueError(f"Overlap detected between called_base_files and og_alignment_file: {', '.join(overlap)}")

    lazy_scaffold = pl.scan_parquet(og_scaffold_file)
    valid_sites = [0, 1, 3, 4]
    new_seqs = {}

    for called_base in called_files:
        sample_id = os.path.basename(called_base).replace("_Called.parquet", "")

        lazy_new = (
            pl.scan_parquet(called_base)
            .filter(pl.col("type").is_in(valid_sites))
            .select(["contig_index", "contig_position", "final_base"])
        )

        lazy_join = (
            lazy_scaffold
            .join(lazy_new, on=['contig_index', 'contig_position'], how="left")
            .with_columns(pl.col("final_base").fill_null("N"))
            .select(["final_base"])
        )

        df = lazy_join.collect()

        new_seqs[sample_id] = "".join(df["final_base"].to_list())


    new_alignment_file = os.path.abspath(new_alignment_file)

    with open(new_alignment_file, "w", encoding="utf-8") as out, open(og_alignment_file, "r", encoding="utf-8") as og:
        for line in og:
            out.write(line)

        for sample_id in called_ids:
            out.write(f">{sample_id}\n")
            for line in textwrap.wrap(new_seqs[sample_id], 80):
                out.write(line + "\n")

if __name__ == "__main__":

    og_alignment_file,og_scaffold_file,called_base_files,new_alignment_file = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]
    append_alignment(og_alignment_file,og_scaffold_file,called_base_files,new_alignment_file)