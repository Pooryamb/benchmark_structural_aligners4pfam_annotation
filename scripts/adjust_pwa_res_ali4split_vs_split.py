import glob
import os
import json
import pandas as pd
from parse_stockholm import parse_stockholm

base_dir = "."
target_stos = glob.glob(f"{base_dir}/data/raw/dbs/pfam_split_target/grp_by_family/PF*.sto")

# Here, the IDs of the target split for each family are extracted
fam_target_ids = {}
for target_path in target_stos:
    sto_dict, _ = parse_stockholm(open(target_path).read())
    fam_id = list(sto_dict.keys())[0]
    target_ids = list(list(sto_dict.values())[0].keys())

    fam_target_ids[fam_id] = target_ids


# Read hmmscan residue alignment information
hmmscan_df = pd.read_csv(f"{base_dir}/data/processed/residue_ali_frac_per_seed_split_vs_split/hmmscan_all.tsv", sep="\t")
hmmscan_df["family"] = hmmscan_df["query"].str.split("-", expand=True)[3]

# Find the tools, residue types, and the families
tools = [os.path.basename(x) for x in glob.glob(f"{base_dir}/tmp/intrafam_residue_alignment_counts/*")]
res_types = [os.path.basename(x) for x in glob.glob(f"{base_dir}/tmp/intrafam_residue_alignment_counts/{tools[0]}/*")]
fams = [os.path.basename(x).replace('.tsv', '') for x in glob.glob(f"{base_dir}/tmp/intrafam_residue_alignment_counts/{tools[0]}/all/*")]

# Find the query target pairs that were used in hmmscan output
query_target_ids_in_hmmscan_output = {}
for family, grp_data in hmmscan_df.groupby("family"):
    # Find the ids of the queries
    query_ids = grp_data["query"]
    # Find the ids of targets
    target_ids = fam_target_ids[family]
    query_target_ids_in_hmmscan_output[family] = (query_ids, target_ids)


for tool in tools:
    for kind in res_types:
        sel_ali_dfs = []
        ali_files = glob.glob(f"{base_dir}/tmp/intrafam_residue_alignment_counts/{tool}/{kind}/*.tsv")
        for ali_path in ali_files:
            fam = os.path.basename(ali_path).replace(".tsv", "")
            all_ali_data = pd.read_csv(ali_path, sep="\t")
            q_t_rep_in_hmmscan = query_target_ids_in_hmmscan_output.get(fam, ([], []))
            desired_ali_data = all_ali_data[(all_ali_data["query"].isin(q_t_rep_in_hmmscan[0])) & (all_ali_data["target"].isin(q_t_rep_in_hmmscan[1]))]
            sel_ali_dfs.append(desired_ali_data)
        sel_ali_rows = pd.concat(sel_ali_dfs)
        sel_ali_rows = sel_ali_rows[sel_ali_rows["all"]!=0] # This is to avoid division by zero
        sel_ali_rows["aligned_fraction"] = sel_ali_rows["correct"]/sel_ali_rows["all"] # Find the percentage of correctly aligned residues

        averaged_ali = sel_ali_rows.groupby("query")["aligned_fraction"].agg("mean").reset_index()  # Calculate the average alignment percentage for a query protein
        averaged_ali.to_csv(f"{base_dir}/data/processed/residue_ali_frac_per_seed_split_vs_split/{tool}_{kind}.tsv", sep="\t", index=None)
