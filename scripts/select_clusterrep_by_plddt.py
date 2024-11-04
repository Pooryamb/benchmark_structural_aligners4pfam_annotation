import pandas as pd

def max_avg_plddt(df):
    max_plddt_afid = df.iloc[df["avg_plddt"].argmax()]['afid']
    df["new_rep"] = max_plddt_afid
    return df

avg_plddt_path = "./data/pfam_avg_plddt.tsv"
clust_path = "./data/dbs/pfam_cif_cut/pfam_clust_cluster.tsv"
output_path = "./data/pfam_clust_reps.tsv"

plddt_df = pd.read_csv(avg_plddt_path, sep="\t", header=None, names=["afid", "avg_plddt"])
clust_df = pd.read_csv(clust_path, sep="\t", header=None, names=["clust_afid", "afid"])
clust_df["clust_afid"] = clust_df["clust_afid"].str.replace(".cif", "")
clust_df["afid"] = clust_df["afid"].str.replace(".cif", "")

clust_plddt = clust_df.merge(plddt_df, on="afid")
clust_new_rep = clust_plddt.groupby("clust_afid").apply(max_avg_plddt, include_groups=False).reset_index(drop=True)
clust_new_rep["new_rep"] = clust_new_rep["new_rep"] + ".cif"
clust_new_rep[["new_rep"]].drop_duplicates().to_csv(output_path, header=None, index=None)

