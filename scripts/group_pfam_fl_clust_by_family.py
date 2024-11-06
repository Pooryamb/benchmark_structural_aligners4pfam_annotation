import pandas as pd
import os

pfam_base_path = "./data/dbs/pfam_fs_cut_clust/"
pfam_fl_base_path ="./data/dbs/pfam_fl_clust/"

df_seeds_h = pd.read_csv(f"{pfam_base_path}/pfam_h.tsv", sep="\t", header=None)
df_seeds_h["db_id"] = df_seeds_h[1].str.split("-", expand=True)[3]
df_seeds_h["fl_id"] = df_seeds_h[1].str.split("-", expand=True)[0]

df_h = pd.read_csv(f"{pfam_fl_base_path}/pfam_fl_h.tsv", sep="\t", header=None, index_col=0)
df_h["unip_id"] = df_h[1].str.split("-", expand=True)[1]

df_seq = pd.read_csv(f"{pfam_fl_base_path}/pfam_fl.tsv", sep="\t", header=None, index_col=0)
df_ss = pd.read_csv(f"{pfam_fl_base_path}/pfam_fl_ss.tsv", sep="\t", header=None, index_col=0)
df_ca = pd.read_csv(f"{pfam_fl_base_path}/pfam_fl_ca.tsv", sep="\t", header=None, index_col=0)

unique_dbids = df_seeds_h["db_id"].unique()


for db_id in unique_dbids:
    fl_ids = df_seeds_h[df_seeds_h["db_id"]==db_id]["fl_id"].unique()
    indices = df_h["unip_id"].isin(fl_ids)

    os.system(f"mkdir -p {pfam_fl_base_path}/grp_by_family/{db_id}")
    output_basename = f"{pfam_fl_base_path}/grp_by_family/{db_id}/{db_id}"

    df_h_sel = df_h.loc[indices, [1]].reset_index(drop=True).reset_index(drop=False)
    df_h_sel.to_csv(f"{output_basename}_h.tsv", sep="\t", header=None, index=None)

    df_seq_sel = df_seq.loc[indices, [1]].reset_index(drop=True).reset_index(drop=False)
    df_seq_sel.to_csv(f"{output_basename}.tsv", sep="\t", header=None, index=None)

    df_ss_sel = df_ss.loc[indices, [1]].reset_index(drop=True).reset_index(drop=False)
    df_ss_sel.to_csv(f"{output_basename}_ss.tsv", sep="\t", header=None, index=None)

    df_ca_sel = df_ca.loc[indices, [1]].reset_index(drop=True).reset_index(drop=False)
    df_ca_sel.to_csv(f"{output_basename}_ca.tsv", sep="\t", header=None, index=None)
