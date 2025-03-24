import json
import pandas as pd

all_seeds = pd.read_csv("./data/raw/dbs/pfam_fs_cut_clust/pfam_h.tsv", sep="\t", header=None, names=["index", "seed_id"]).drop(columns=["index"])
all_seeds["size"] = all_seeds["seed_id"].str.split("-", expand=True)[2].astype(int) - all_seeds["seed_id"].str.split("-", expand=True)[1].astype(int) + 1

all_res_dict = {}
for index, row in all_seeds.iterrows():
    all_res_dict[row["seed_id"]] = list(range(1, row["size"]+1))
with open("./tmp/residue_features/all.json", 'w') as ofile:
    json.dump(all_res_dict, ofile, indent=4)
