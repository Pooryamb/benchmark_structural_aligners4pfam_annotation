import pickle
import glob
import os
import multiprocessing as mp
import pandas as pd
import numpy as np


script_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(script_dir, "..")

with open(f"{base_dir}/data/processed/parsed_pssm.pkl", 'rb') as ifile:
    pssm_data = pickle.load(ifile)


k_df = pd.DataFrame(pssm_data["k_values"].items(), columns=["pfam", "k"])
lamda_df = pd.DataFrame(pssm_data["lamda_value"].items(), columns=["pfam", "lamda"])
del pssm_data

def calc_bitscore(input_path):
    "Reads the alignments scored by psiblast scheme and then calculates the bitscore representative"
    ali_df = pd.read_csv(input_path, sep="\t")
    ali_df = ali_df.merge(k_df, left_on="t_family", right_on="pfam").drop(columns=["pfam"])
    ali_df = ali_df.merge(lamda_df, left_on="t_family", right_on="pfam").drop(columns=["pfam"])
    ali_df["bitscore_rep"] = ali_df["lamda"] * ali_df["pssm_raw_score"] - np.log2(ali_df["k"]) # As the lamda value is identical for all families, we simply ignore it
    ali_df = ali_df.drop(columns=["lamda", "k"])
    ali_df.to_csv(input_path.replace("_rescored.tsv", "_bitscore.tsv"), sep="\t", index=None)

input_paths = glob.glob(f"{base_dir}/tmp/alis/split_pf_seq/*_rescored.tsv")
with mp.Pool(4) as pool: #processes=mp.cpu_count()
    batch_data_list = pool.map(calc_bitscore, input_paths)
