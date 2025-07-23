import os
import sys
import glob
import pickle
import pandas as pd
import numpy as np
import numba as nb
from fasta2dict import fasta2dict
from map_residues_numba import map_residues_numba

script_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(script_dir, "..")

# Find the mapping between each seed and its family consensus
from fasta2dict import fasta2dict
fasta_paths = glob.glob(f"{base_dir}/data/raw/dbs/pfam_split_target/grp_by_family/*.fasta")

mapping2cons = {}
for fasta_path in fasta_paths:
    parsed_fasta = fasta2dict(fasta_path)
    for seq_id, seq in parsed_fasta.items():
        if seq_id.endswith("_consensus"):
            cons_seq = seq
        else:
            res_mapping = pd.DataFrame(map_residues_numba(seq.encode(), cons_seq.encode()), columns=["seed_coord", "cons_coord"])
            mapping2cons[seq_id] = res_mapping

with open(f"{base_dir}/data/processed/seed2cons_mapping.pkl", 'wb') as ofile:
    pickle.dump(mapping2cons, ofile)
