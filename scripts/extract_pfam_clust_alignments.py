import gzip
import pandas as pd
from parse_stockholm import parse_stockholm, write_parsed_sto
data_path = "./data/"



# Step 1: read the sequences from PfamSDB. For both PfamDB and PfamSDB, sequences are read in a dictionary format
# The key of each dictionary is the family accession, and its values are dictionaries whose keys are fasta headers
# And values are sequences
pfams_fam_seq_dict = {}
with open(f"{data_path}/raw/dbs/pfam_fs_cut_clust/pfam_h.tsv") as h_file, \
     open(f"{data_path}/raw/dbs/pfam_fs_cut_clust/pfam.tsv") as seq_file:
    while True:
        header_line = h_file.readline()
        if not(header_line.strip()):
            break
        header = header_line.strip().split()[1]
        family = header.split("-")[3]
        seq = seq_file.readline().strip().split()[1]
        fam_seqs = pfams_fam_seq_dict.get(family, {})
        fam_seqs[header] = seq
        pfams_fam_seq_dict[family] = fam_seqs


# Step 2: Read sequences from PfamDB.
pfam_fam_seq_dict, acc2id = parse_stockholm(gzip.open(f'{data_path}/raw/Pfam-A.seed.gz', 'rt').read())

# Align PfamSDB sequences with each other based on Pfam alignments
aligned_pfams_seqs = {}
for fam, fam_seqs_pfams in pfams_fam_seq_dict.items():
    fam_seqs_pfam = pfam_fam_seq_dict[fam]
    for pfams_id, pfams_seq in fam_seqs_pfams.items():
        for pfam_id, pfam_seq in fam_seqs_pfam.items():
            if pfam_seq.replace(".", '').upper() == pfams_seq:
                pfams_aligned_dict = aligned_pfams_seqs.get(fam, {})
                pfams_aligned_dict[pfams_id] = pfam_seq
                aligned_pfams_seqs[fam] = pfams_aligned_dict


write_parsed_sto(aligned_pfams_seqs, acc2id, f"{data_path}/processed/aligned_pfam_clust.sto")
