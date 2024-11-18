import argparse
import os
import pandas as pd


parser = argparse.ArgumentParser(description="""This script selects a random subset of Pfam seeds. 
The output of this script will be used to make a foldseek subdb""")

parser.add_argument("--fraction", type=float, help="The fraction of seeds to be sampled. This must be a number between 0 and 1")
parser.add_argument("--number", type=int, help="The numbere of seeds to be sampled")

args = parser.parse_args()

if not (args.number or args.fraction):
    raise Exception("Please either provide fraction or number")


src_df = pd.read_csv("./data/pfam_clust_reps.tsv", sep="\t", header=None)
if args.fraction:
    sampled_df = src_df.sample(frac=args.fraction, random_state=42)
else:
    sampled_df = src_df.sample(n=args.number, random_state=42)

sampled_df.to_csv("./data/pfam_sample_reps.tsv", sep="\t", header=None, index=None)
