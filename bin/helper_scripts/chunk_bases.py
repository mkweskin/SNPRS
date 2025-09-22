import sys
import polars as pl
import os
import json
import traceback
import psutil
import glob
import shutil

def chunk_sample(scaffold_parquet, sample_parquet, sample_id,chunk_json, temp_directory):
    try:
        
        with open(chunk_json, "r") as f:
            raw_coords = json.load(f)            
        chunk_coords = [(int(i), (int(rng[0]), int(rng[1]))) for i, rng in raw_coords]  
        
        frames = []
        for cid, (start, end) in chunk_coords:
            df = pl.DataFrame({"row_idx": range(start, end)})
            df = df.with_columns(pl.lit(cid).alias("chunk_id"))
            frames.append(df)

        chunk_map = pl.concat(frames, rechunk=True)
        
        scaffold_df = (
            pl.scan_parquet(scaffold_parquet)
            .select(["contig_index", "contig_position"])
            .with_row_index("row_idx")
            .with_columns(pl.col("row_idx").cast(pl.Int64))
        )
        
        scaffold_df = scaffold_df.join(chunk_map.lazy(), on="row_idx", how="left")     
           
        sample_df = (
            pl.scan_parquet(sample_parquet)
            .filter(pl.col("type").is_in([0, 1, 4, 5]))
            .select([
                "contig_index",
                "contig_position",
                pl.when(pl.col("type") == 5)
                .then(pl.lit("N"))
                .otherwise(pl.col("final_base"))
                .alias("final_base"),
            ])
        )

        joined = (
            scaffold_df
            .join(sample_df, on=["contig_index", "contig_position"], how="left")
            .select(["chunk_id", "final_base"])
        ).collect()
        

        chunk_files = []
        for cid, (start, end) in chunk_coords:
            
            chunk_file = os.path.join(
                temp_directory,
                os.path.basename(sample_parquet).replace(".parquet", f"_Chunk_{cid}.parquet")
            )
            chunk_files.append(chunk_file)

            (
                joined
                .filter(pl.col("chunk_id") == cid)
                .select("final_base")
                .write_parquet(chunk_file, compression="snappy")
            )

        output = {"chunk_files": chunk_files}
        print(json.dumps(output))
    
    except Exception as e:
        print(f"❌ Error: {e}")
        traceback.print_exc()

if __name__ == "__main__":
    scaffold_parquet, sample_parquet, sample_id,chunk_json, temp_directory = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    chunk_sample(scaffold_parquet, sample_parquet, sample_id,chunk_json, temp_directory)
