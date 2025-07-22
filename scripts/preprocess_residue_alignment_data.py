import json
import glob
import os
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import pickle
import pandas as pd
import numpy as np
from sklearn.utils import resample
from scipy.interpolate import interp1d
from joblib import Parallel, delayed


script_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(script_dir, "../")

def cumulative_aligned_ratio(df):
    """Given a dataframe containing the fraction of correctly aligned residues for each query, it returns a cumulative
    dataframe showing the fraction of queries whose average aligned fraction is above a specific threshold"""
    sorted_df = df.sort_values(by="aligned_fraction", ascending=False)
    sorted_df["numb"] = 1
    total_size = sorted_df.shape[0]
    sorted_df["cum_frac"] = sorted_df["numb"].cumsum() / total_size
    return sorted_df[["aligned_fraction", "cum_frac"]]

tool_name = {"fs3di": "Foldseek (3Di)", "fs": "Foldseek", "mm": "MMseqs", "rs": "Reseek", "tm": "TM-align", "hmmscan": "Hmmscan"}

def calc_auc_once(df):
    """Receiving the cumulative alignment ratio dataframe, it will output the area under the curve"""
    y_axis = np.insert(df["aligned_fraction"], 0, 1)
    x_axis = np.insert(df["cum_frac"], 0, 0)
    auc = np.trapz( y_axis, x_axis) 
    return auc


def to_parallel_auc(i, df):
    """This was defined to only do the time-consuming part with one parallelized function"""
    sampled_df = resample(df, replace=True, random_state=i)
    cum_samp_df = cumulative_aligned_ratio(sampled_df)
    return calc_auc_once(cum_samp_df)


def preprocess_res_ali_plot(df, n_bootstrap=1000, confidence=0.95):
    """Receiving the dataframe containing the alignment fraction for each query, it will output the y, x, auc, and confidence intervals for auc"""
    with ProcessPoolExecutor() as executor:
        bs_aucs = list(executor.map(partial(to_parallel_auc, df=df), range(n_bootstrap)))

#    bs_aucs = []
#    for i in range(n_bootstrap):
#        bs_aucs.append(to_parallel_auc(i, df))
    
    lower_range = (1 - confidence)/2
    upper_range = 1 - lower_range
    
    lower_ci_auc = np.quantile(bs_aucs, lower_range ) # Stores the low-bound CI for auc
    upper_ci_auc = np.quantile(bs_aucs, upper_range ) # Stores the high-bound CI for auc
    
    cum_df = cumulative_aligned_ratio(df) # First, find the cumulative plot showing the fraction of queries whose fraction of correctly aligned residues is above a threshold
    overall_auc = calc_auc_once(cum_df)
    y_axis = np.insert(cum_df["aligned_fraction"], 0, 1) # 1 is inserted in the beginning to make sure AUC calculation considers y=1 for small x values
    x_axis = np.insert(cum_df["cum_frac"], 0, 0)
    return {"x_axis":x_axis, "y_axis":y_axis, "auc": overall_auc, "lower_ci_auc": lower_ci_auc, "upper_ci_auc": upper_ci_auc} 


# From now on, the residue alignment efficiency for each kind of residue is calculated

## Overall residue alignment block
paths = sorted(glob.glob(f"{base_dir}/data/processed/residue_ali_frac_per_seed_split_vs_split/*_all.tsv"))
overall_res_ali_dict = {}
for path in paths:
    tool_abbreviation = os.path.basename(path).split("_")[0]
    df = pd.read_csv(path, sep="\t")
    plot_data = preprocess_res_ali_plot(df)
    overall_res_ali_dict[tool_abbreviation] = plot_data


## Conserved residue alignment block
paths = sorted(glob.glob(f"{base_dir}/data/processed/residue_ali_frac_per_seed_split_vs_split/*_conserved.tsv"))
conserved_res_ali_dict = {}
for path in paths:
    tool_abbreviation = os.path.basename(path).split("_")[0]
    df = pd.read_csv(path, sep="\t")
    plot_data = preprocess_res_ali_plot(df)
    conserved_res_ali_dict[tool_abbreviation] = plot_data


## Active site alignment block
paths = sorted(glob.glob(f"{base_dir}/data/processed/residue_ali_frac_per_seed_split_vs_split/*_active_site.tsv"))
active_site_res_ali_dict = {}
for path in paths:
    tool_abbreviation = os.path.basename(path).split("_")[0]
    df = pd.read_csv(path, sep="\t")
    plot_data = preprocess_res_ali_plot(df)
    active_site_res_ali_dict[tool_abbreviation] = plot_data
    

## Binding site alignment block
paths = sorted(glob.glob(f"{base_dir}/data/processed/residue_ali_frac_per_seed_split_vs_split/*_binding_site.tsv"))
binding_site_res_ali_dict = {}
for path in paths:
    tool_abbreviation = os.path.basename(path).split("_")[0]
    df = pd.read_csv(path, sep="\t")
    plot_data = preprocess_res_ali_plot(df)
    binding_site_res_ali_dict[tool_abbreviation] = plot_data
    

## Low confidence conserved_sites block

seeds_with_low_confidence_conserved_sites = pd.read_csv(f"{base_dir}/data/processed/low_confidence_conserved_residues.tsv", sep="\t")

paths = sorted(glob.glob(f"{base_dir}/data/processed/residue_ali_frac_per_seed_split_vs_split/*_conserved.tsv"))
lccs_res_ali_dict = {}
for path in paths:
    tool_abbreviation = os.path.basename(path).split("_")[0]
    df = pd.read_csv(path, sep="\t")
    df = df[df["query"].isin(seeds_with_low_confidence_conserved_sites["seed_id"])]
    plot_data = preprocess_res_ali_plot(df)
    lccs_res_ali_dict[tool_abbreviation] = plot_data

all_residue_alignment_data = {"low_confidence_conserved_sites": lccs_res_ali_dict,
                              "overall": overall_res_ali_dict, 
                              "conserved": conserved_res_ali_dict, 
                              "active_sites": active_site_res_ali_dict, 
                              "binding_sites": binding_site_res_ali_dict
                              }

with open(f"{base_dir}/data/processed/residue_alignment_performance.pkl", 'wb') as f:
    pickle.dump(all_residue_alignment_data, f)