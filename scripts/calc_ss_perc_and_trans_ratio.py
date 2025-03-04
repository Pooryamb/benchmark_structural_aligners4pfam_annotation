import pandas as pd
from fasta2dict import fasta2dict

data_path = "./data/"

ss_fasta_path = f"{data_path}/processed/ss_pfam.fasta"
seeds_tsv_path = f"{data_path}/raw/dbs/pfam_fs_cut_clust/pfam_h.tsv"
fasta_dict = fasta2dict(ss_fasta_path)
seeds_df = pd.read_csv(seeds_tsv_path, sep="\t", header=None)
seeds_df[["unip_id", "start", "end", "family"]] = seeds_df[1].str.split("-", expand=True)
seeds_df["start"] = seeds_df["start"].astype(int)
seeds_df["end"] = seeds_df["end"].astype(int)

def count_transitions(secondary_structure_sequence):  
    transitions = 0  
    for i in range(1, len(secondary_structure_sequence)):  
        if secondary_structure_sequence[i] != secondary_structure_sequence[i - 1]:  
            transitions += 1  
    return transitions  

ss_fracs = []
for index, row in seeds_df.iterrows():
    fl_ss = fasta_dict[row["unip_id"]]
    dom_ss = fl_ss[row["start"]-1: row["end"]]
    dom_len = row["end"] - row["start"] + 1
    h_frac = dom_ss.count("H") / dom_len
    e_frac = dom_ss.count("E") / dom_len
    c_frac = dom_ss.count("C") / dom_len
    tr_count = count_transitions(dom_ss)
    norm_tr = tr_count / dom_len
    ss_fracs.append([row[1], h_frac, e_frac, c_frac, norm_tr])
    
ss_frac_df = pd.DataFrame(ss_fracs, columns=["seed_id", "h_frac", "e_frac", "c_frac", "len_norm_tr_count"])

ss_frac_df.to_csv(f"{data_path}/processed/ss_info_pfam.tsv", sep="\t", index=None)
