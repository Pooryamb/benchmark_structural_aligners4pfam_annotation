import glob
import pickle
import os
import pandas as pd
from read_text_blocks import read_text_blocks
from parse_stockholm import parse_stockholm

script_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(script_dir, "..")

# First, we find the consensus sequence for the built MSA

hmm_path = f"{base_dir}/data/raw/dbs/pfam_split_target/pfam.hmm"

con_seq_dict = {}
for hmm_block in read_text_blocks(hmm_path):
    # parse the Pfam accession
    pfam_acc = [x.split()[1] for x in hmm_block if x.startswith("ACC")][0]
    # find the line where hmm starts
    hmm_start_ind = [x for x in range(len(hmm_block)) if hmm_block[x].startswith("HMM          A")][0]
    # Now, we only select the lines containing the MSA and HMM positions
    pos_lines = hmm_block[hmm_start_ind+5:-1:3]
    consensus_seq = "".join([x.split()[-4] for x in pos_lines])
    consensus_seq = consensus_seq.upper()
    con_seq_dict[pfam_acc] = consensus_seq


# Now, we load the mapping between MSA and HMM columns
with open(f"{base_dir}/data/processed/hmm2msa_mapping.pkl", 'rb') as ifile:
    mapping_data = pickle.load(ifile)


def fmt_con_based_on_mapping(con, mapping_df):
    """This function takes a consensus sequence and also a dataframe containing the mapping between the positions of an
    MSA and the HMM. Then, it will return the consensus sequence formatted based on the mapping df"""
    last_pos = 0
    aas = []
    for i, aa in enumerate(con):
        pos_in_aligned_con = mapping_df.loc[i]["msa_col"]
        gaps_num = pos_in_aligned_con - len(aas) - 1
        aas = aas + gaps_num * ["-"]
        aas.append(aa)
    return "".join(aas)


# The consensus sequence is aligned with the MSA of its family
aligned_cons_seq_dict = {}
for key, fam_mapping_df in mapping_data.items():
    aligned_seq = fmt_con_based_on_mapping(con_seq_dict[key], mapping_data[key])
    aligned_cons_seq_dict[key] = aligned_seq

# Make the MSA including the Consensus sequence from Pfam MSA

for family in mapping_data:
    aligned_cons_seq = aligned_cons_seq_dict[family]
    fasta_in_dict = {f"{family}_concensus": aligned_cons_seq}
    sto_path = f"{base_dir}/data/raw/dbs/pfam_split_target/grp_by_family/{family}.sto"
    parsed_sto = parse_stockholm(open(sto_path).read())[0]
    seq_dict = list(parsed_sto.values())[0]
    fasta_in_dict = {x: y.replace(".", "-") for x,y in seq_dict.items()}
    fasta_out_path = sto_path.replace(".sto", ".fasta")
    with open(fasta_out_path, 'w') as ofile:
        ofile.write(f">{family}_consensus\n{aligned_cons_seq}\n")
        cons_len = len(aligned_cons_seq)  # The length of each seed must be equal to the length of the consensus
        for key, value in fasta_in_dict.items():
            ofile.write(f">{key}\n{value[:cons_len]}\n")
