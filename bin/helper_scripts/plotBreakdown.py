import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Plot SNPRS data")
    
    parser.add_argument("--in", dest="snp_file", type=str, required=True,help="Path to SNP data")
    parser.add_argument("--out", dest="output_file", type=str, required=True,help="Path to output file")
    
    return parser.parse_args()

args = parse_args()
snp_file = os.path.abspath(args.snp_file)
out_file = os.path.abspath(args.output_file)

valid_exts = [".pdf", ".png"]
ext = os.path.splitext(out_file)[1].lower()

if ext not in valid_exts:
    print(f"Error: Output file must end with .pdf or .png. Got '{out_file}'.")
    sys.exit(1)

df_results = pd.read_csv(snp_file, sep="\t")

species_order = sorted(df_results["Species"].unique())

df_results["Species"] = pd.Categorical(df_results["Species"], categories=species_order, ordered=True)

def plot_all_samples(df_results,save_path):
    samples = df_results["Sample_ID"].unique()

    n = len(samples)
    ncols = 5
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4 * nrows), squeeze=False)

    for ax, sample_id in zip(axes.flatten(), samples):
        df = df_results[df_results["Sample_ID"] == sample_id].copy()

        df = df.sort_values("Proportion", ascending=True)

        ax.errorbar(
            df["Proportion"],
            df["Species"],
            xerr=[df["Proportion"] - df["CI_lower"],
                  df["CI_upper"] - df["Proportion"]],
            fmt='o',
            capsize=4,
            linewidth=1.4
        )
        ax.set_title(sample_id)
        ax.set_xlabel("Proportion")
        ax.set_ylabel("Species")

    # Turn off empty subplots
    for i in range(len(samples), nrows * ncols):
        axes.flatten()[i].axis("off")

    plt.tight_layout()
    plt.savefig(save_path)

plot_all_samples(df_results, save_path=out_file)


