import glob
import os
import pickle
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
from sklearn.utils import resample
from scipy.interpolate import interp1d
from headers import headers

script_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(script_dir, "../data")
base_dir = os.path.join(script_dir, "../")
                        
def add_label_col(df, is_hmmscan=False):
    """Adds two columns named q_pfam, and t_pfam, and also a column showing if they are the same for a row"""
    df["q_pfam"] = df["query"].str.split("-", expand=True)[3]
    if is_hmmscan:
        df["t_pfam"] = df["t_acc"]
    else:
        df["t_pfam"] = df["target"].str.split("-", expand=True)[3]
    df["pfam_label"] =  df["t_pfam"]==df["q_pfam"]
    return df

def add_evalue_bin(df, is_tm=False):
    "Adds a new column for binning the e-value of the hits, for TM-align, -100* length normalized tm-score is used instead"
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    if is_tm:
        df["evalue_bin"] = -np.round(100 * df['evalue'])
    else:
        df["evalue_bin"] = np.where(df['evalue'] == 0, -1000, np.round(np.log2(df['evalue'])))
    return df

# Find the precision and recall at each e-value threshold:
q_id_df = pd.read_csv(f"{data_dir}/raw/dbs/pfam_split_query/pfam_h.tsv", sep="\t", header=None)
total_q_num = q_id_df.shape[0]  # The total number of queries for recall calculation

def calc_precision_recall_vs_evalue(df):
    """Sorts the alignment table based on its e-value and then at each e-value threshold, it calculates the precistion and recall"""
    df = df.sort_values(by="evalue_bin")
    df["cum_tp"] = df["pfam_label"].cumsum()
    df = df.reset_index(drop=True)
    df = df.reset_index(names="row_num")
    df["row_num"] = df["row_num"] + 1
    df["cum_fp"] = df["row_num"] - df["cum_tp"]
    df["precision"] = df["cum_tp"]/df["row_num"]
    df["recall"] = df["cum_tp"]/total_q_num
    df["f1"] = 2 * df["precision"] * df["recall"] / (df["precision"] + df["recall"])
    return df.drop_duplicates("evalue_bin", keep="last")

def calc_maxf1_precision_recall_evalue(df):
    """Given a dataframe containing the F1 score, it returns the maximum F1 along with precision, recall, and the evalue leading to such a F1"""
    max_f1 = df["f1"].max()
    max_f1_arg = df["f1"].argmax()
    max_f1_row = df.iloc[max_f1_arg]
    return  { "max_f1": max_f1, 
              "max_f1_precision": max_f1_row["precision"], 
              "max_f1_recall": max_f1_row["recall"],
              "max_f1_evalue_bin": max_f1_row["evalue_bin"]}

def _bootstrap_f1_once(df, i):
    """Bootstraps the dataframe and calculates the maximum f1 score. For the maximum
    f1-score, it reports the precision, recall, and the evalue_bin"""
    bs_df = resample(df, replace=True, random_state=i)
    bs_df_f1 = calc_precision_recall_vs_evalue(bs_df)
    return calc_maxf1_precision_recall_evalue(bs_df_f1)["max_f1"]
    

def bootstrap_f1(df, n_bootstrap=1000, confidence=0.95):
    """It finds the CI for the maximum F1 score"""
    
    """
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(_bootstrap_f1_once, df, i)
            for i in range(n_bootstrap)
        ]
        bootstrap_f1_values = [f.result() for f in futures] 
    """
    
    bootstrap_f1_values = []
    
    for i in range(n_bootstrap):
        bootstrap_f1_values.append(_bootstrap_f1_once(df, i)) 
    # Calculate confidence intervals
    lower_range = (1 - confidence)/2
    upper_range = 1 - lower_range
    
    lower_f1 = np.quantile(bootstrap_f1_values, lower_range ) # Stores the low-bound CI for F1
    upper_f1 = np.quantile(bootstrap_f1_values, upper_range ) # Stores the high-bound CI for F1

    return {"f1_ci_lower": lower_f1, "f1_ci_upper": upper_f1}


# Copy all the paths related to each tool into a dictionary

tsv_files = sorted(glob.glob(f"{base_dir}/tmp/first_hits/*.tsv"))
method_paths_dict = {}
for path in tsv_files:
    method = os.path.basename(path).split("_B")[0]
    method_paths = method_paths_dict.get(method, [])
    method_paths.append(path)
    method_paths_dict[method] = method_paths

method_df_dict = {method: pd.concat([pd.read_csv(x, sep="\t") for x in method_paths_dict[method]]) for method in method_paths_dict} # Reads every thing in dataframes
method_df_dict = {x:add_label_col(y, "hmmscan" in x) for x,y in method_df_dict.items()}  # Adds label at family can clan level
method_df_dict = {x:add_evalue_bin(y, "tm" in x) for x,y in method_df_dict.items()}      # Adds e-value bins


f1_ci_dict = {x:bootstrap_f1(y) for x,y in method_df_dict.items()} # Find the F1 CIs

prec_recall = {x: calc_precision_recall_vs_evalue(y) for x,y in method_df_dict.items()}  # Calculates precision and recall for each tool
f1_dict = {x: calc_maxf1_precision_recall_evalue(y) for x, y in prec_recall.items()}

trim_first_rows_frac = 0.05
prec_recall = {x: y[y["row_num"] >= trim_first_rows_frac * total_q_num] for x,y in prec_recall.items()}  # Get rid of the first 5% because it is noisy

all_f1_data = {x: {**f1_ci_dict[x], **f1_dict[x], "precision_vs_recall": prec_recall[x]} for x in method_df_dict.keys()}
with open(f"{data_dir}/processed/f1_data.pkl", 'wb') as file:
    pickle.dump(all_f1_data, file)