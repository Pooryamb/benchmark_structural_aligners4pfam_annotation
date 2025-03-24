import pandas as pd

def max_avg_plddt(df):
    max_plddt_afid = df.iloc[df["avg_plddt"].argmax()]['seed_id']
    df["new_rep"] = max_plddt_afid
    return df

avg_plddt_path = "./data/processed/pfam_avg_plddt.tsv"
clust_path = "tmp/alis/pfam_clust/pfam_clust_cluster.tsv"
#./data/dbs/pfam_cif_cut/pfam_clust_cluster.tsv"
output_path = "./data/processed/pfam_clust_reps.tsv"

plddt_df = pd.read_csv(avg_plddt_path, sep="\t")
clust_df = pd.read_csv(clust_path, sep="\t", header=None, names=["clust_seed_id", "seed_id"])
clust_df["clust_seed_id"] = clust_df["clust_seed_id"].str.replace(".cif", "")
clust_df["seed_id"] = clust_df["seed_id"].str.replace(".cif", "")

clust_plddt = clust_df.merge(plddt_df, on="seed_id")
clust_new_rep = clust_plddt.groupby("clust_seed_id").apply(max_avg_plddt, include_groups=False).reset_index(drop=True)
clust_new_rep["len"] = clust_new_rep["seed_id"].str.split("-", expand=True)[2].astype(int) - clust_new_rep["seed_id"].str.split("-", expand=True)[1].astype(int)
clust_new_rep = clust_new_rep[clust_new_rep["len"] >= 10]
clust_new_rep["new_rep"] = clust_new_rep["new_rep"] + ".cif"
clust_new_rep[["new_rep"]].drop_duplicates().to_csv(output_path, header=None, index=None)

