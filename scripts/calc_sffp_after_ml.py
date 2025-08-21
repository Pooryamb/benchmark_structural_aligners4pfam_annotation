import glob
import sys
import os
import multiprocessing as mp
import pickle
import pandas as pd
from calc_sffp4split_alis_rescored import *


if __name__== "__main__":
    foldseek_paths = [f"{base_dir}/tmp/alis/split_pf_seq/fs_pref_B{i}_ml2.tsv" for i in range(3,17)]     ##### Single batch used for training
    reseek_paths =   [f"{base_dir}/tmp/alis/split_pf_seq/reseek_sens_B{i}_ml2.tsv" for i in range(3,17)] ##### Single batch used for training
    cols2read = ["query", "target", "evalue", "bitscore_rep", "ml_score"]
    tps_bef_first_fp = {
                        "foldseek_ml": find_tps_count_bef_first_fp_for_filepaths(foldseek_paths,cols2read, col2sortby="ml_score", sort_ascending=False),
                        "foldseek": find_tps_count_bef_first_fp_for_filepaths(foldseek_paths,cols2read, col2sortby="evalue", sort_ascending=True),
                        "reseek_ml": find_tps_count_bef_first_fp_for_filepaths(reseek_paths, cols2read, col2sortby="ml_score", sort_ascending=False),
                        "reseek": find_tps_count_bef_first_fp_for_filepaths(reseek_paths, cols2read,col2sortby="evalue", sort_ascending=True)
                        }

#,
#
#   "foldseek_bitscore": find_tps_count_bef_first_fp_for_filepaths(foldseek_paths, cols2read, col2sortby="bitscore_rep", sort_ascending=False),
#   ,
#   "reseek": find_tps_count_bef_first_fp_for_filepaths(reseek_paths, cols2read,col2sortby="evalue", sort_ascending=True),
#   "reseek_bitscore": find_tps_count_bef_first_fp_for_filepaths(reseek_paths, cols2read, col2sortby="bitscore_rep", sort_ascending=False)
    sffp_dir_path = f"{base_dir}/tmp/sffp_rescoring"
    os.makedirs(sffp_dir_path, exist_ok=True)
    for key, value in tps_bef_first_fp.items():
        value.to_csv(f"{sffp_dir_path}/{key}_exc_b12.tsv", sep="\t", index=None)  ##### Exc b1 or b12
