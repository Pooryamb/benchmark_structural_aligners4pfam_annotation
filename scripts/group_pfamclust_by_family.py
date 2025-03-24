import pandas as pd
import os


def read_fasta_to_df(fasta_file_path, exclude_cif_ext=True):
    """Reads a fasta file in a dataframe format where the first column stores the header of sequences and
    the second column stores its sequence. If the inclusion_list is not empty, it will only read the
    sequences whose header is in that inclusion_list"""
    rows = []
    with open(fasta_file_path) as file:
        for line in file:
            if not(line.strip()):
                break
            if line.startswith(">"):
                header = line.strip().lstrip(">")
                if exclude_cif_ext:
                    header = header.replace(".cif", "")
            else:
                seq = line.strip()
                rows.append([header, seq])

    return pd.DataFrame(rows)


data_dir = "./data/"
seq_df = read_fasta_to_df(f"{data_dir}/raw/dbs/pfam_cif_cut_clust/pfam.fasta")
ss_df = read_fasta_to_df(f"{data_dir}/raw/dbs/pfam_cif_cut_clust/pfam_ss.fasta")[[1]].rename(columns={1:2})
ca_df = read_fasta_to_df(f"{data_dir}/raw/dbs/pfam_cif_cut_clust/pfam_ca.fasta")[[1]].rename(columns={1:3})

concat_df = pd.concat([seq_df, ss_df, ca_df], axis=1)
del seq_df, ss_df, ca_df
concat_df["family"] = concat_df[0].str.split("-", expand=True)[3]

for fam, fam_df in concat_df.groupby("family"):
    out_dir = f"{data_dir}/raw/dbs/pfam_cif_cut_clust/grp_by_family/{fam}"
    os.makedirs(out_dir, exist_ok=True)
    fam_df[[0]].reset_index(drop=True).reset_index().to_csv(f"{out_dir}/{fam}_h.tsv", sep="\t", index=None, header=None)
    fam_df[[1]].reset_index(drop=True).reset_index().to_csv(f"{out_dir}/{fam}.tsv", sep="\t", index=None, header=None)
    fam_df[[2]].reset_index(drop=True).reset_index().to_csv(f"{out_dir}/{fam}_ss.tsv", sep="\t", index=None, header=None)
    fam_df[[3]].reset_index(drop=True).reset_index().to_csv(f"{out_dir}/{fam}_ca.tsv", sep="\t", index=None, header=None)
