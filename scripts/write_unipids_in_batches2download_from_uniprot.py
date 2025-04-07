import pandas as pd
import math

prots = pd.read_csv("./data/raw/dbs/pfam_fs_cut_clust/pfam_h.tsv", sep="\t", header=None)[[1]]

unip_accessions = prots[1].str.split("-", expand=True)[0].unique()

batch_size = 90_000 # UniProt accepts  up to 100_000 queries at a time
batch_num = math.ceil(len(unip_accessions) / 90_000)
for i in range(batch_num):
    selected_ids = unip_accessions[i*batch_size: (i+1)*batch_size]
    with open(f"./tmp/extra/unip_ids_b{i}.txt", 'w') as ofile:
        ofile.write("\n".join(selected_ids))
