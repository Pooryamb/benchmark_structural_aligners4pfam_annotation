import os
import glob
from concurrent.futures import ProcessPoolExecutor
import pickle
import numpy as np
import pandas as pd
from sklearn.utils import resample
from scipy.interpolate import interp1d
from joblib import Parallel, delayed
from sens_up2_first_fp import get_cum_sens_from_1st_occ_df, get_tps_frac_bef_1st_fp_from_1st_occ_df, cum_sens_dist, auc_tps_bef_1st_fp


fig_dir = "./figures"
os.makedirs(fig_dir, exist_ok=True)
### This part reads all tables from a tool and make a single dataframe out of it.

data_dir = "./data"
first_occ_dict = {}

all_tsv_files = glob.glob(f"{data_dir}/processed/first_label_occ/*_B*.tsv")
tools = set(list([os.path.basename(x).split("_B")[0] for x in all_tsv_files])) # The tools will be reseek, fs_cut, cif_cut, mm, and tm
file_paths_tool_dict = {x: [y for y in all_tsv_files if x in y] for x in tools} 

# As the first occurences were found in parallel computations, the following lines merge the result for each batch
for tool, tool_tsv_paths in file_paths_tool_dict.items():
    first_occ_list = []
    for tsv_path in tool_tsv_paths:
        tsv_df = pd.read_csv(tsv_path, sep="\t")
        first_occ_list.append(tsv_df)
    first_occ_dict[tool] = pd.concat(first_occ_list)

### Now, we make two dataframes, one for storing the cumulative sffp plot (cum_sens_up2first_fp_dict_fam), 
### and one for storing the fraction of TPs before first FP for each query (frac_sens_up2first_fp_dict_fam)
### A function is defined to make such dictionaries in parallel mode.

def calc_cumsffp_and_frac(tool):
    """For the tool, it will calculate the cumulative SFFP plot for both family and clan level, it will also find the fraction
    of TPs before the first FP for each tool"""
    cum_sens_fam = get_cum_sens_from_1st_occ_df(first_occ_dict[tool], "tp_bef_fp_frac_pfam")
    cum_sens_clan = get_cum_sens_from_1st_occ_df(first_occ_dict[tool], "tp_bef_fp_frac_clan")
    frac_sens = get_tps_frac_bef_1st_fp_from_1st_occ_df(first_occ_dict[tool])
    return tool, cum_sens_fam, cum_sens_clan, frac_sens


# Find the fraction of TPs before first FP, and also the cumulative SFFP plot for each tool
cum_sens_up2first_fp_dict_fam = {}
cum_sens_up2first_fp_dict_clan = {}
frac_sens_up2first_fp_dict_fam = {}


with ProcessPoolExecutor() as executor:
    results = executor.map(calc_cumsffp_and_frac, first_occ_dict.keys())
    for tool, cum_sens_fam, cum_sens_clan, frac_sens in results:
        cum_sens_up2first_fp_dict_fam[tool] = cum_sens_fam
        cum_sens_up2first_fp_dict_clan[tool] = cum_sens_clan
        frac_sens_up2first_fp_dict_fam[tool] = frac_sens

### The functions defined for CI estimations by bootstrapping

def auc_1st_fp_func(df):
    return auc_tps_bef_1st_fp(df, "tp_bef_fp_frac_pfam")

def _bootstrap_single_stat(i, df, x_values, target_col="tp_bef_fp_frac_pfam"):
    """It first bootstraps the data and then calculates the AUC for sffp along with estimated sffp for each x value. Here, x is the "fraction of queries". 
    it takes an i as the random state, df, which is a dataframe showing the the sffp for each query, x_values, were explained earlier
    target_col will be the column used for the whole analysis. It can be
    "tp_bef_fp_frac_pfam" or "tp_bef_fp_frac_clan" """
    boot_sample = df.sample(n=len(df), replace=True, random_state=i)
    cum_df = cum_sens_dist(boot_sample)
    # It never has repetetive numbers cum_df = cum_df.drop_duplicates(["cum_sum_frac"])
    x_axis, y_axis = cum_df["cum_sum_frac"], cum_df[target_col]
    interp_sffp = interp1d(x_axis, y_axis,
            bounds_error=False, fill_value="extrapolate", assume_sorted=False
        )
    interp_y = interp_sffp(x_values)
    auc = np.trapz(y_axis, x_axis)
    return auc, interp_y

def bootstrap_statistic(df, x_values, n_bootstrap=1000, confidence=0.95):
    """It will bootstrap for SFFP_AUC and SFFP plot calculations"""
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(_bootstrap_single_stat, i, df, x_values)
            for i in range(n_bootstrap)
        ]
        bootstrap_auc = [f.result()[0] for f in futures]
        bootstrap_sffp_data = pd.DataFrame([f.result()[1] for f in futures])

    # Calculate confidence intervals
    lower_range = (1 - confidence)/2
    upper_range = 1 - lower_range
    
    sffp_lower_ci = bootstrap_sffp_data.quantile(lower_range)   # Stores the low-bound CI for the SFFP plot
    sffp_upper_ci = bootstrap_sffp_data.quantile(upper_range)   # Stores the high-bound CI for the SFFP plot
    
    lower_auc = np.percentile(bootstrap_auc, lower_range * 100) # Stores the low-bound CI for the AUC of SFFP plot
    upper_auc = np.percentile(bootstrap_auc, upper_range * 100) # Stores the high-bound CI for the AUC of SFFP plot

    return {
        'auc_ci_lower': lower_auc,
        'auc_ci_upper': upper_auc,
        'sffp_ci_lower': sffp_lower_ci,
        'sffp_ci_upper': sffp_upper_ci
    }

## Store the CI information for all tools in a dictionary

ci_fam_dict = {}

for tool in file_paths_tool_dict:
    x_values = cum_sens_up2first_fp_dict_fam[tool]["cum_sum_frac"]
    ci_fam_dict[tool] = bootstrap_statistic(frac_sens_up2first_fp_dict_fam[tool], x_values)


all_sffp_data = {
                 "sffp": cum_sens_up2first_fp_dict_fam, 
                 "sffp_clan": cum_sens_up2first_fp_dict_clan,
                 "auc_sffp": {tool: auc_1st_fp_func(frac_sens_up2first_fp_dict_fam[tool]) for tool in file_paths_tool_dict},
                 "ci_info": ci_fam_dict,
                 "frac_info": frac_sens_up2first_fp_dict_fam
                }


with open(f'{data_dir}/processed/sffp.pkl', 'wb') as f:
    pickle.dump(all_sffp_data, f)
