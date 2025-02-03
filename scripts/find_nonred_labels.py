import argparse
import os
import pandas as pd


parser = argparse.ArgumentParser(description="This script takes "
    "the output alignment of a tool as input and outputs the sensitivity "
    "up to the first FP at both Pfam and Clan level. It also requires "
    "a file showing the associations between Pfam ID and Clan ID.")

parser.add_argument("--input", type=str, required=True, help="The path to the input alignment file")
parser.add_argument("--output", type=str, help="The first occurence of each label "
"will be saved here. It shows the first Pfam or Clan level label")
parser.add_argument("--pfam_clan_info", default="./data/raw/pfam_clan_info.tsv", type=str, help= "The path to the file relating pfam to clan")
parser.add_argument("--remove_cif_ext", type=str, default="Auto", help= "Do you want the cif extension? True/False/Auto")

args = parser.parse_args()

if not args.output:
    args.output = "data/processed/first_label_occ/" + os.path.basename(args.input)

clan_info_path = args.pfam_clan_info
ipr_clan_df = pd.read_csv(clan_info_path, sep="\t")

if args.remove_cif_ext == "True":
    remove_cif = True
elif args.remove_cif_ext == "False":
    remove_cif = False
elif args.remove_cif_ext == "Auto":
    if ".cif" in open(args.input).readline().split()[0]:
        remove_cif = True
    else:
        remove_cif = False

chunksize = 10_000
selected = []
i = 0
for chunk in pd.read_csv(args.input, sep="\t", header=None, chunksize=10_000):
    if remove_cif:
        chunk[0] =  chunk[0].str.replace(".cif", "")
        chunk[1] =  chunk[1].str.replace(".cif", "")
    chunk["q_pfam"] = chunk[0].str.split("-", expand=True)[3]
    chunk["t_pfam"] = chunk[1].str.split("-", expand=True)[3]
    chunk = chunk.reset_index(names=['row_num'])

    chunk = chunk.merge(ipr_clan_df, left_on="q_pfam", right_on="pfam").rename(columns={"clan": "q_clan"}).drop(columns=["pfam"])
    chunk = chunk.merge(ipr_clan_df, left_on="t_pfam", right_on="pfam").rename(columns={"clan": "t_clan"}).drop(columns=["pfam"])
    chunk["pfam_label"] = (chunk["q_pfam"] == chunk["t_pfam"])
    chunk["clan_label"] = (chunk["q_clan"] == chunk["t_clan"])
    i += 1
    first_occ = chunk.drop_duplicates([0, "pfam_label", "clan_label"])
    selected.append(first_occ)


selected_df = pd.concat(selected)
selected_df = selected_df.drop_duplicates([0, "pfam_label", "clan_label"])

selected_df.to_csv(args.output, sep="\t", header=None, index=False)
