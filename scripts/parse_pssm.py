import os
import sys
import glob
import pickle
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import numba as nb

script_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(script_dir, "..")

def parse_pssm(pssm_path):
    """receiving a pssm path, it will parse it and return the pssm along with k and lamda"""
    df = pd.read_csv(pssm_path, delim_whitespace=True, header=None, engine="python", skiprows=3, skipfooter=3)
    pssm_df = df.iloc[:, 2:22].rename(columns = {x:x-2 for x in df.columns}) # Read the tabular part and rename the columns
    k, lamda = open(pssm_path).read().strip().strip().split()[-2:]
    k = float(k)
    lamda = float(lamda)
    return pssm_df, k, lamda

# Parse the pssm files

pssm_dict = {}
k_dict = {}
lamda_dict = {}
pssm_paths = glob.glob(f"{base_dir}/data/raw/dbs/pfam_split_target/grp_by_family/*.pssm")
for pssm_path in pssm_paths:
    family = os.path.basename(pssm_path).replace(".pssm", "")
    pssm_dict[family], k_dict[family], lamda_dict[family] = parse_pssm(pssm_path)

# Find the mapping between each amino acid and the its header number in pandas

with open(pssm_paths[0]) as file:
    aas = file.readlines()[2].strip().split()[:20]
aa_inds_dict = {aa: i for i,aa in enumerate(aas)}

all_pssm_data = {"pssm_profiles": pssm_dict, "k_values": k_dict, "lamda_value": lamda_dict, "aa_ind": aa_inds_dict}
with open(f"{base_dir}/data/processed/parsed_pssm.pkl", 'wb') as ofile:
    pickle.dump(all_pssm_data, ofile)
