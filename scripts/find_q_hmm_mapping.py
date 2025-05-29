import glob
import pickle
import pandas as pd
from read_text_blocks import read_text_blocks
from find_pw_correspondence import find_pw_correspondence


base_dir = "."
large_files = glob.glob(f"{base_dir}/tmp/alis/fam_alis/hmmscan/PF*_large_file.txt")


def parse_ali_block_in_hmm(block):
    """Takes an expanded alignment block as input and returns the summarized pairwise alignment as output. The output
    includes the query accession, hmm_accession, q_start, q_end, t_start, t_end, query sequence, and target sequence."""
    # In case there are multiple alignments between a query and an HMM, one with the minimum e-value is selected
    block_parts = "\n".join(block).split("conditional E-value: ")[1:]
    eval2blocks = {float(x.split()[0]): "\n".join(x.split("\n")[1:-1]) for x in block_parts}
    min_e_value = min(eval2blocks.keys())
    min_e_ind = list(eval2blocks.keys()).index(min_e_value)
    min_e_block = list(eval2blocks.values())[min_e_ind].split("\n")
    seq_lines = min_e_block[0:-1]
    q_lines = [x for x in seq_lines[2::5] if x.strip()]
    t_lines = [x for x in seq_lines[0::5] if x.strip()]
    q_seq = "".join([x.split()[2] for x in q_lines])
    t_seq = "".join([x.split()[2] for x in t_lines])
    q_id = q_lines[0].split()[0]
    t_id = t_lines[0].split()[0]
    q_start = int(q_lines[0].split()[1])
    q_end = int(q_lines[-1].split()[-1])
    t_start = int(t_lines[0].split()[1])
    t_end = int(t_lines[-1].split()[-1])
    return q_id, t_id, q_start, q_end, t_start, t_end, q_seq, t_seq


def get_seq_hmm_mapping(parsed_block):
    "This takes the parsed block as input and returns the mapping between query sequence and hmm in a dataframe"
    q_id, t_id, q_start, q_end, t_start, t_end, q_seq, t_seq = parsed_block
    mappings = find_pw_correspondence(q_seq, t_seq, q_start, t_start)
    return pd.DataFrame(mappings, columns=["q_col", "hmm_col"])


block_start = "  Alignments for each domain:"
block_end = "Internal pipeline statistics summary:"
q_hmm_mappings = {}
for large_file_path in large_files:

    for block in read_text_blocks(large_file_path, block_start=block_start, block_end=block_end):
        parsed_block = parse_ali_block_in_hmm(block)
        q_hmm_mappings[parsed_block[0]] = get_seq_hmm_mapping(parsed_block)

output_path = f"{base_dir}/data/processed/q2hmm_mapping.pkl"
with open(output_path, 'wb') as file:
    pickle.dump(q_hmm_mappings, file)
