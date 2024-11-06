import pandas as pd
import os

base_path = "./data/dbs/pfam_fs_cut_clust/"
df_h = pd.read_csv(f"{base_path}/pfam_h.tsv", sep="\t", header=None, index_col=0)
df_seq = pd.read_csv(f"{base_path}/pfam.tsv", sep="\t", header=None, index_col=0)
df_ss = pd.read_csv(f"{base_path}/pfam_ss.tsv", sep="\t", header=None, index_col=0)
df_ca = pd.read_csv(f"{base_path}/pfam_ca.tsv", sep="\t", header=None, index_col=0)

df_h["db_id"] = df_h[1].str.split("-", expand=True)[3]
for grp in df_h.groupby("db_id"):
    grp_index = grp[1].index
    grp_id = grp[0]
    os.system(f"mkdir -p {base_path}/grp_by_family/{grp_id}")
    output_basename = f"{base_path}/grp_by_family/{grp_id}/{grp_id}"

    df_h_sel = df_h.loc[grp_index, [1]].reset_index(drop=True).reset_index(drop=False)
    df_h_sel.to_csv(f"{output_basename}_h.tsv", sep="\t", header=None, index=None)

    df_seq_sel = df_seq.loc[grp_index, [1]].reset_index(drop=True).reset_index(drop=False)
    df_seq_sel.to_csv(f"{output_basename}.tsv", sep="\t", header=None, index=None)

    df_ss_sel = df_ss.loc[grp_index, [1]].reset_index(drop=True).reset_index(drop=False)
    df_ss_sel.to_csv(f"{output_basename}_ss.tsv", sep="\t", header=None, index=None)

    df_ca_sel = df_ca.loc[grp_index, [1]].reset_index(drop=True).reset_index(drop=False)
    df_ca_sel.to_csv(f"{output_basename}_ca.tsv", sep="\t", header=None, index=None)
