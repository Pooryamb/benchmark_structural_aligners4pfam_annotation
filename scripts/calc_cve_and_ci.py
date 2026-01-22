import glob
import os
import pickle
import pandas as pd
import numpy as np
from sklearn.utils import resample
from scipy.interpolate import interp1d
from multiprocessing import Pool, cpu_count
from cve import cve

script_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(script_dir, "../data")


def bootstrap_once(i, q_grouped_df, cov_x, q_count):
    """Selects a random sample of queries and interpolates the fpepq based on the cov_x values fed into it
    q_count is needed because the number of queries in bootstraped samples will not be equal to the unique number of queries.
    As for sampling the dataframe, we need to group the dataframe, we group it only once so that in each sampling, we just
    combine the groups that we need.
    """

    selected_qs = resample(pd.Series(q_grouped_df.keys()), replace=True, random_state=i) # Select a sample of the queries
    query_counts = selected_qs.value_counts()   # Find the number of times each query has been sampled

    sampled_df_lst = []
    for q, count in query_counts.items():
        sampled_df_lst = sampled_df_lst + count * [q_grouped_df[q]]    # Copy the data for each sampled query as many as it should be selected
    boot_sample = pd.concat(sampled_df_lst)   # Make a dataframe of all selected data
    
    cve_boot_df = cve(boot_sample, q_count=q_count).drop_duplicates(["cov_pfam"]).drop_duplicates(["fpepq_pfam"]).reset_index(drop=True) # Find the CVE plot for the sampled data. As the number of queries is important and cannot be found from the sampled data, it has to be provided to the CVE function separately.

    # The CVE is estimated on the x values used in the orignial plot
    interp_cve = interp1d(
        cve_boot_df["cov_pfam"], cve_boot_df["fpepq_pfam"],
        bounds_error=False, fill_value="extrapolate", assume_sorted=False
    )
    
    boot_interp_fpepq = interp_cve(cov_x)

    return boot_interp_fpepq



def bootstrap_cve_parallel(df, cov_x, n_bootstrap=1000, n_jobs=None):
    """This runs the function in parallel"""
    q_count = df["query"].nunique()
    q_grouped_df = {grp_name: grp_df for grp_name, grp_df in df.groupby("query")}

    with Pool(processes=cpu_count()) as pool:
        results = pool.starmap(
            bootstrap_once,
            [(i, q_grouped_df, cov_x, q_count) for i in range(n_bootstrap)]
        )

    fpepq_results = pd.DataFrame([r for r in results if r is not None])

    return fpepq_results


# First, read all data from e-value bins and concatenate them to make a single dataframe
evalue_bins_path = f"{data_dir}/processed/evalue_bins/"
file_paths = glob.glob(f"{evalue_bins_path}/*.tsv")
tools = sorted(list(set(os.path.basename(x).split("_B")[0] for x in file_paths)))
tools_paths = {x: sorted(glob.glob(f"{evalue_bins_path}/{x}_B*.tsv")) for x in tools}
# Files of each batch are read and concatenated together.
tools_df = {}
for tool in tools:
    dfs = [pd.read_csv(x, sep="\t") for x in tools_paths[tool]]
    tools_df[tool] = pd.concat(dfs)


upper_fpepq_threshold = 10 # We want Sensitivity at FPEPQ up to 10
cve_dict = {}
# Find the upper and lower quantiles
ci_value = 0.95
lower_range = (1 - ci_value)/2
upper_range = 1 - lower_range

for tool, df in tools_df.items():
    print(f"processing {tool}")
    cve_df = cve(df).drop_duplicates(["cov_pfam"]).drop_duplicates(["fpepq_pfam"]).reset_index(drop=True) # Calculate the CVE for the orignial data

    interpolated_rel_inv = interp1d(cve_df["fpepq_pfam"], cve_df["cov_pfam"], bounds_error=False, fill_value="extrapolate") # Train a model to find the coverage corresponding to FPEPQ=10
    
    # Add 10 to the FPEPQ values and select FPEPQ values up to 10
    fpepq_values = pd.concat([cve_df["fpepq_pfam"], pd.Series([upper_fpepq_threshold])]).drop_duplicates().sort_values().reset_index(drop=True) 
    fpepq_values = fpepq_values[fpepq_values<=upper_fpepq_threshold]
    
    cov_x = interpolated_rel_inv(fpepq_values) # To find the cov_x whose FPEPQ is equal to 10

    boot_strapped_covs  = bootstrap_cve_parallel(df, cov_x) 

    lower_ci_cov = boot_strapped_covs.quantile(lower_range)
    upper_ci_cov = boot_strapped_covs.quantile(upper_range)

    cve_dict[tool] = [cov_x, fpepq_values, lower_ci_cov, upper_ci_cov]

with open(f'{data_dir}/processed/cve_ci.pkl', 'wb') as f:
    pickle.dump(cve_dict, f)
