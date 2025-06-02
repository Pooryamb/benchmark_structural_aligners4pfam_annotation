import os
import pandas as pd
import numpy as np
from fasta2dict import fasta2dict
from dict2fasta import dict2fasta
from parse_stockholm import parse_stockholm

base_dir = "."

whole_pfam_clust = fasta2dict(f"{base_dir}/data/raw/dbs/pfam_cif_cut_clust/pfam.fasta")
whole_pfam_clust = {x.replace(".cif", ""):y for x,y in whole_pfam_clust.items()}
t_db = pd.read_csv(f"{base_dir}/data/raw/dbs/pfam_split_target/pfam_h.tsv", sep="\t", header=None)
all_q_ids = set(whole_pfam_clust.keys()) - set(t_db[1])
query_fasta = {x:y for x,y in whole_pfam_clust.items() if x in all_q_ids}

grouped_fasta = {}

for seed_id, seq in query_fasta.items():
    fam = seed_id.split("-")[-1]
    seq_dict = grouped_fasta.get(fam, {})
    seq_dict[seed_id] = seq
    grouped_fasta[fam] = seq_dict

out_dir = f"{base_dir}/data/raw/dbs/pfam_split_query/grp_by_family/"
os.makedirs(out_dir, exist_ok=True)

for fam, seq_dict in grouped_fasta.items():
    dict2fasta(seq_dict, f"{out_dir}/{fam}.fasta")
