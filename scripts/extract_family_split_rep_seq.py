import pickle
from read_text_blocks import read_text_blocks



# This script finds the mapping between hmm and its source MSA

base_dir = "."
hmm_path = f"{base_dir}/data/raw/dbs/pfam_split_target/pfam.hmm"


all_seqs = {}
for hmm_block in read_text_blocks(hmm_path):
    # parse the Pfam accession
    pfam_acc = [x.split()[1] for x in hmm_block if x.startswith("ACC")][0]
    # find the line where hmm starts
    hmm_start_ind = [x for x in range(len(hmm_block)) if hmm_block[x].startswith("HMM          A")][0]
    # Now, we only select the lines containing the representative sequence states
    pos_lines = hmm_block[hmm_start_ind+5:-1:3]
    rep_seq = "".join([x.split()[-4] for x in pos_lines]).upper()
    all_seqs[pfam_acc] = rep_seq


with open(f"{base_dir}/data/processed/pfam_fam_rep_seq.fasta", 'w') as ofile:
    for fam_id, seq in all_seqs.items():
        ofile.write(f">{fam_id}\n{seq}\n")

