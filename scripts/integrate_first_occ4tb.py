import argparse
import os
import glob
import pandas as pd

parser = argparse.ArgumentParser(description="This script combines the split files containing the row_number of the first occurence of each label, to have a single output file.")
parser.add_argument("--input_dir", type=str, required=True, help="The path to the directory containing the split files. It is assumed that all files are sorted alphabetically.")
parser.add_argument("--output", type=str, help="The directory where the split files will be stored. If not provided, the input directory will be used, a filename with the same name as the input directory will be created (+ .tsv)")
parser.add_argument("--batch_size", type=int, default=10_000_000, help="The maximum number of lines per split file (default: 10 million)")

args = parser.parse_args()

input_dir = args.input_dir
input_files = sorted(glob.glob(os.path.join(input_dir, "*.tsv")))
output = args.output if args.output else input_dir + ".tsv"

row_num_corrected_dfs = []
for idx, input_file in enumerate(input_files):
    df = pd.read_csv(input_file, sep='\t')
    df["row_num"] += idx * args.batch_size
    row_num_corrected_dfs.append(df)

combined_df = pd.concat(row_num_corrected_dfs, ignore_index=True)
combined_df = combined_df.drop_duplicates(["query", "pfam_label", "clan_label"], keep="first")

combined_df.to_csv(output, sep='\t', index=False)
