import pickle
import pandas as pd
from parse_stockholm import parse_stockholm


base_dir = "."
pfam_fam_seq_dict, acc2id = parse_stockholm(open(f"{base_dir}/data/processed/aligned_pfam_clust.sto").read())
seed2msa_mapping = {}

for pf, pf_alis in pfam_fam_seq_dict.items():
    ali_seq_len = len(list(pf_alis.values())[0])
    for seed_id, seed_seq in pf_alis.items():
        mappings = []
        seed_cursor = 0
        for msa_col_1 in range(ali_seq_len):
            if seed_seq[msa_col_1].isalnum():
                seed_cursor += 1
                mappings.append([msa_col_1 + 1, seed_cursor])

        seed2msa_mapping[seed_id] = pd.DataFrame(mappings, columns=["msa_col", "seed_col"])

with open(f"{base_dir}/data/processed/seed2msa_mapping.pkl", 'wb') as ofile:
    pickle.dump(seed2msa_mapping, ofile)
