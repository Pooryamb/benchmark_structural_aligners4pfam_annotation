import pandas as pd


seeds_df = pd.read_csv("./data/raw/dbs/pfam_fs_cut_clust/pfam_h.tsv", sep="\t", header=None)
unip_ids = seeds_df[1].str.split("-", expand=True)[0].unique()
with open("./tmp/extra/urls.txt", 'w') as ofile:
    for acc in unip_ids:
        ofile.write(f"gs://public-datasets-deepmind-alphafold-v4/AF-{acc}-F1-model_v4.cif\n")

