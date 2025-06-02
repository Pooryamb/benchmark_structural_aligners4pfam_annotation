import pickle
import pandas as pd
from read_text_blocks import read_text_blocks



# This script finds the mapping between hmm and its source MSA

base_dir = "."
hmm_path = f"{base_dir}/data/raw/dbs/pfam_split_target/pfam.hmm"


all_mappings = {}
for hmm_block in read_text_blocks(hmm_path):
    # parse the Pfam accession
    pfam_acc = [x.split()[1] for x in hmm_block if x.startswith("ACC")][0]
    # find the line where hmm starts
    hmm_start_ind = [x for x in range(len(hmm_block)) if hmm_block[x].startswith("HMM          A")][0]
    # Now, we only select the lines containing the MSA and HMM positions
    pos_lines = hmm_block[hmm_start_ind+5:-1:3]
    pos_mapping = [[x.split()[0], x.split()[-5]] for x in pos_lines]
    mapping_df = pd.DataFrame(pos_mapping, columns=["hmm_col", "msa_col"])
    mapping_df = mapping_df[["hmm_col", "msa_col"]].astype(int)
    #mapping_df["pf"] = pfam_acc
    all_mappings[pfam_acc] = mapping_df
    #all_mappings.append(mapping_df)
#pfam_hmm2msa_mapping = pd.concat(all_mappings)


#pfam_hmm2msa_mapping.to_csv(f"{base_dir}/data/processed/hmm2msa_mapping.tsv", sep="\t", index=None)
with open(f"{base_dir}/data/processed/hmm2msa_mapping.pkl", 'wb') as ofile:
    pickle.dump(all_mappings, ofile)
