import argparse
import os
import pandas as pd
from headers import headers
from label_pf_clan_chunk import label_pf_clan_chunk

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

selected = []
for chunk in label_pf_clan_chunk(args.input):
    first_occ = chunk.drop_duplicates(["query", "pfam_label", "clan_label"])
    selected.append(first_occ)

selected_df = pd.concat(selected)
selected_df = selected_df.drop_duplicates(["query", "pfam_label", "clan_label"])

selected_df.to_csv(args.output, sep="\t", index=False)
