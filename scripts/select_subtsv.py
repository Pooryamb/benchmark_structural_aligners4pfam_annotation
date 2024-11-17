import os
import argparse
import pandas as pd


parser = argparse.ArgumentParser(description="""This script as input receives the path to a tsv
file, and also a file containing some IDs, it also has an option for removing .cif extensions. 
Then, it makes new database that is a subset of the original one""")

parser.add_argument("--input_basename", type=str, required=True, help="The path to the tsv file. It should not contain the .tsv extension")
parser.add_argument("--list2sel", type=str, required=True, help="The path to the file containing the list of the files to select.")
parser.add_argument("--remove_cif_extension", type=str, help="Whether the .cif extension should be removed from the members of the list")
parser.add_argument("--output_basename", type=str, required=True, help="The path to the output tsv file basenames")


args = parser.parse_args()

os.system(f"mkdir -p  {os.path.dirname(args.output_basename)}")

tosel_df = pd.read_csv(f"{args.list2sel}", sep="\t", header=None)
if args.remove_cif_extension == "True":
    tosel_df[0] = tosel_df[0].str.replace(".cif", "")

tsv_df_exts = {"seq": ".tsv", "ss": "_ss.tsv", "ca": "_ca.tsv", "h": "_h.tsv"}

h_df = pd.read_csv(f"{args.input_basename}_h.tsv", sep="\t", header=None, index_col=0)
indices = h_df[h_df[1].isin(tosel_df[0])].index


for db_type, ext in tsv_df_exts.items():
    df = pd.read_csv(f"{args.input_basename}{ext}", sep="\t", header=None, index_col=0)
    df_sel = df.loc[indices, 1]
    df_sel = df_sel.reset_index(drop=True).reset_index(drop=False)
    df_sel.to_csv(f"{args.output_basename}{ext}", sep="\t", header=None, index=None)
