import pickle
import glob
import os
import json
import pandas as pd

base_dir = "."
hmm2msa_dict = {}
special_residues_files = glob.glob(f"{base_dir}/tmp/residue_features/*.json")
special_residues_dict = {os.path.basename(x).replace(".json", ""): json.loads(open(x).read()) for x in special_residues_files}

#load msa information
with open(f'{base_dir}/data/processed/hmm2msa_mapping.pkl', 'rb') as hmm2msa_handle, \
     open(f'{base_dir}/data/processed/q2hmm_mapping.pkl', 'rb') as q2hmm_handle, \
     open(f'{base_dir}/data/processed/seed2msa_mapping.pkl', 'rb') as seed2msa_handle:
    hmm2msa = pickle.load(hmm2msa_handle)
    q2hmm = pickle.load(q2hmm_handle)
    seed2msa = pickle.load(seed2msa_handle)

all_file_path = f'{base_dir}/data/processed/residue_ali_frac_per_seed_split_vs_split/hmmscan_all.tsv'
special_files_paths = {kind: f"{base_dir}/data/processed/residue_ali_frac_per_seed_split_vs_split/hmmscan_{kind}.tsv" for kind in special_residues_dict.keys()}

with open(all_file_path, 'w') as all_sum_handler:
    special_file_handlers = {kind: open(file_path, 'w') for kind, file_path in special_files_paths.items()}
    header_line = "query\taligned_fraction\n"
    all_sum_handler.write(header_line)
    for f in special_file_handlers.values():
        f.write(header_line)

    for q_name, q_prof_mapping in q2hmm.items():
        q_fam = q_name.split("-")[-1]
        prof2msa = hmm2msa[q_fam]
        q_msa_hmmer = q_prof_mapping.merge(prof2msa, on="hmm_col")[["q_col", "msa_col"]]
        q_msa_gs = seed2msa[q_name] # gs stands for the gold standard
        q_len = int(q_name.split("-")[2]) - int(q_name.split("-")[1]) + 1 #The query length
        # Check overal residue alignment
        gs_hmmer_shared = q_msa_hmmer.merge(q_msa_gs, left_on=["q_col", "msa_col"], right_on=["seed_col", "msa_col"])
        shared_mappings_size = gs_hmmer_shared.shape[0]
        all_sum_handler.write(f"{q_name}\t{shared_mappings_size/q_len}\n")
        # Now check residue alignment for special residues
        for res_kind, res_list in special_residues_dict.items():
            seed_special_residues = res_list.get(q_name, [])
            if seed_special_residues:
                correctly_aligned_special_residues = gs_hmmer_shared[gs_hmmer_shared["q_col"].isin(seed_special_residues)]
                special_file_handlers[res_kind].write(f"{q_name}\t{len(correctly_aligned_special_residues)/len(seed_special_residues)}\n")

    for f in special_file_handlers.values():
        f.close()
