import gzip
import pandas as pd
from parse_stockholm import parse_stockholm, write_parsed_sto
data_path = "./data/"


#A function for finding the sequence identity between two aliged sequences. It is assumed that alignments
# are in the format available in Pfam database.
def pident_of_aligned_seqs(seq1, seq2):
    identical_cols = sum([1 for i in range(len(seq1)) if (seq1[i].upper()==seq2[i].upper() and (seq1[i]!="." or seq2[i]!="."))])
    total_length = sum([1 for i in range(len(seq1)) if (seq1[i]!="." or seq2[i]!=".")])
    return identical_cols/total_length * 100

aligned_pfams_seqs, acc2id = parse_stockholm(open(f"{data_path}/processed/aligned_pfam_clust.sto").read())

seq_ids = []
for family, fam_dict in aligned_pfams_seqs.items():
    for seq_id_1, seq_1 in fam_dict.items():
        for seq_id_2, seq_2 in fam_dict.items():
            seq_ids.append([seq_id_1, seq_id_2, pident_of_aligned_seqs(seq_1, seq_2)])

seq_ids_df = pd.DataFrame(seq_ids)
seq_ids_df = seq_ids_df[seq_ids_df[0]!=seq_ids_df[1]]

summary_df = seq_ids_df.groupby(0)[2].agg("mean").reset_index().rename(columns={0:"seed_id", 2: "avg_intra_fam_pident"})
summary_df.to_csv(f"{data_path}/processed/avg_intra_fam_pident.tsv", sep="\t", index=None)
