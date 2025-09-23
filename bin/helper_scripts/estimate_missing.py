import polars as pl
import json
import sys
import os
import numpy as np

def find_elbow_index(missing, counts):
    p1 = np.array([missing[0], counts[0]])
    p2 = np.array([missing[-1], counts[-1]])

    line_vec = p2 - p1
    line_vec_norm = line_vec / np.sqrt(np.sum(line_vec**2))

    distances = []
    for i in range(len(missing)):
        p = np.array([missing[i], counts[i]])
        vec = p - p1
        proj_len = np.dot(vec, line_vec_norm)
        proj_point = p1 + proj_len * line_vec_norm
        dist = np.sqrt(np.sum((p - proj_point) ** 2))
        distances.append(dist)

    return int(np.argmax(distances))

def estimate_missing(site_file,pass_site_type_file,temp_directory):
    
    sites = pl.scan_parquet(site_file).with_row_index("row_nr")
    
    with open(pass_site_type_file, "r", encoding="utf-8") as f:
        pass_site_type_rows = set(json.load(f))
    
    mask = pl.col("row_nr").is_in(list(pass_site_type_rows))
    missing_values = sites.filter(mask).collect()["missing"].to_numpy()
    
    if len(missing_values) == 0:
        return 0
    
    max_obs = int(np.max(missing_values))

    missing_thresholds = np.arange(0, max_obs + 1)
    counts = np.array([np.sum(missing_values <= t) for t in missing_thresholds])

    elbow_idx = find_elbow_index(missing_thresholds, counts)
    ideal_missing = missing_thresholds[elbow_idx]

    print(str(ideal_missing))

if __name__ == "__main__":

    site_file,pass_site_type_file,temp_directory = sys.argv[1],sys.argv[2],sys.argv[3]
    estimate_missing(site_file,pass_site_type_file,temp_directory)
