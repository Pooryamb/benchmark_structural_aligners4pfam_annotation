import glob
import os
import pandas as pd
import numpy as np
from fasta2dict import fasta2dict

qdb_size = 20000
batch_num = 16
base_dir = "."
db_path = f"{base_dir}/data/raw/dbs/"

src_db_basename = f"{db_path}/pfam_cif_cut_clust/pfam"
def remove_cif_ext(seq_dict):
    return {x.replace(".cif", ""): y for x,y in seq_dict.items()}

# Read the fasta file of each database file. You must have converted the dataframes into the fasta file in advance.
p_seq = fasta2dict(f"{src_db_basename}.fasta")
p_seq = remove_cif_ext(p_seq)
ss_seq = fasta2dict(f"{src_db_basename}_ss.fasta")
ss_seq = remove_cif_ext(ss_seq)
ca_seq = fasta2dict(f"{src_db_basename}_ca.fasta")
ca_seq = remove_cif_ext(ca_seq)

seeds = pd.DataFrame({"seed_id" : list(p_seq.keys())})
seeds["family"] = seeds["seed_id"].str.split("-", expand=True)[3]
seeds = seeds.sample(frac=1, random_state=42) # Randomize the hits

# Members of each family will be splitted into two halves. The first half will be less comprehensive than the second half
seeds_first_half = []
seeds_second_half = []

for grp_name, grp_df in seeds.groupby("family"):
    first_half_row_num = grp_df.shape[0] // 2
    second_half_row_num = grp_df.shape[0] - first_half_row_num
    seeds_first_half.append(grp_df.head(first_half_row_num))
    seeds_second_half.append(grp_df.head(second_half_row_num))

seeds_first_half = pd.concat(seeds_first_half).reset_index(drop=True)
seeds_query_df = seeds_first_half.sample(n=qdb_size, random_state=42)
seeds_target_df = pd.concat(seeds_second_half).reset_index(drop=True)

# The accessions of each half are stored for profile making purposes
query_db_members_path = f"{base_dir}/tmp/extra/pfam_split_query_members.txt"
target_db_members_path = f"{base_dir}/tmp/extra/pfam_split_target_members.txt"
seeds_query_df.to_csv(query_db_members_path, sep="\t", index=None)
seeds_target_df.to_csv(target_db_members_path, sep="\t", index=None)

def write_tsv(acc_df, out_path):
    """Takes a dataframe containing the seed_ids as input and writes 4 tsv files in the out_path"""
    with open(f"{out_path}.tsv", 'w') as seq_file, \
         open(f"{out_path}_h.tsv", 'w') as h_file, \
         open(f"{out_path}_ss.tsv", 'w') as ss_file, \
         open(f"{out_path}_ca.tsv", 'w') as ca_file:
        i = 0
        for index, row in acc_df.iterrows():
            _ = h_file.write(f"{i}\t{row.seed_id}\n")
            _ = seq_file.write(f"{i}\t{p_seq[row.seed_id]}\n")
            _ = ss_file.write(f"{i}\t{ss_seq[row.seed_id]}\n")
            _ = ca_file.write(f"{i}\t{ca_seq[row.seed_id]}\n")
            i += 1

# Write the target database
tdb_path = f"{base_dir}/data/raw/dbs/pfam_split_target/pfam"
write_tsv(seeds_target_df, tdb_path)

# Write the query databases

write_tsv(seeds_query_df, f"{base_dir}/data/raw/dbs/pfam_split_query/pfam")
splitted_q = np.array_split(seeds_query_df, batch_num)
for i in range(batch_num):
    subdf = splitted_q[i]
    qdb_path = f"{base_dir}/data/raw/dbs/pfam_split_query/B{i+1}/pfam"
    write_tsv(subdf, qdb_path)
