import pandas as pd
from parse_stockholm import parse_stockholm, write_parsed_sto

base_dir = "."

target_pfam_members = pd.read_csv(f"{base_dir}/tmp/extra/pfam_split_target_members.txt", sep="\t") # Read the file containing the members of the target database
parsed_sto, acc_mapping = parse_stockholm(open(f"{base_dir}/data/processed/aligned_pfam_clust.sto").read()) # Read the stockholm file of all members of the clustered pfam
# Only select the seed members that were part of the target
target_db_dict = {}
for family_id, fam_df in target_pfam_members.groupby("family"):
    family_dict = {}
    for index, row in fam_df.iterrows():
        seed_id = row["seed_id"]
        if seed_id in parsed_sto[family_id]:
            family_dict[seed_id] = parsed_sto[family_id][seed_id]
    target_db_dict[family_id] = family_dict

def rm_gaps_from_sto(parsed_sto):
    """Please note that by far, the sequence alignments have been extracted from the Pfam MSAs. Now, we want to
    remove the columns which have all gaps in the selected sub-database"""
    edited_sto = {}
    for fam, fam_dict in parsed_sto.items():
        seq_aas = {seed_id: [] for seed_id in fam_dict.keys()}
        msa_len = len(list(fam_dict.values())[0])
        for i in range(msa_len):
            ith_aa = {seed_seq[i] for seed_seq in fam_dict.values()}
            if ith_aa == {"."}:
                continue
            for seed_id in fam_dict.keys():
                seq_aas[seed_id].append(fam_dict[seed_id][i])
        edited_sto[fam] = {seed_id: "".join(seq_aas[seed_id]) for seed_id in fam_dict.keys()}
    return edited_sto

# Remove columns that are only gaps
gap_less_sto = rm_gaps_from_sto(target_db_dict)
write_parsed_sto(gap_less_sto, acc_mapping, f"{base_dir}/data/raw/dbs/pfam_split_target/pfam.sto")
