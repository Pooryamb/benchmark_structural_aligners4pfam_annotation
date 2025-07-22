import json
import glob
import os
import pickle
import pandas as pd
import numpy as np
from sklearn.utils import resample
from scipy.interpolate import interp1d
import xlsxwriter
from sklearn.model_selection import train_test_split  
from sklearn.ensemble import RandomForestRegressor  
from sklearn.metrics import mean_squared_error  


script_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(script_dir, "../")
data_dir = f"{base_dir}/data/"


def train_model(data_df):
    """Trains a random forest for each model and shows the importance of the features"""
    # Assuming the last 8 columns are features and the second column is the query  
    X = data_df.iloc[:, 1:-1]     # Features (starting from second column till one to the last columns)  
    y = data_df.iloc[:, -1]       # query (second column)
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    # Initialize the model  
    rf_model = RandomForestRegressor(n_estimators=100, random_state=42)  
    
    # Fit the model  
    rf_model.fit(X_train, y_train)
    # Make predictions  
    y_pred_train = rf_model.predict(X_train)
    y_pred_test = rf_model.predict(X_test)  
    
    # Calculate Mean Squared Error  
    mse_train = mean_squared_error(y_train, y_pred_train)
    mse_test = mean_squared_error(y_test, y_pred_test)
    
    importances = rf_model.feature_importances_  
    # Create a DataFrame for visualization  
    importance_df = pd.DataFrame({'Feature': X.columns, 'Importance': importances})  
    importance_df = importance_df.sort_values(by='Importance', ascending=False)
    
    return {"model":rf_model, "mse_train": mse_train, "mse_test": mse_test, "feature_importance": importance_df}

# Find seed characteristics

data_dir = f"{base_dir}/data/"
pi_df = pd.read_csv(f"{data_dir}/processed/avg_intra_fam_pident.tsv", sep="\t") #pi means percentage identity
ss_info_df = pd.read_csv(f"{data_dir}/processed/ss_info_pfam.tsv", sep="\t")
cn_df = pd.read_csv(f"{data_dir}/processed/avg_contact_num.tsv", sep="\t")
plddt_df = pd.read_csv(f"{data_dir}/processed/pfam_avg_plddt.tsv", sep="\t")
plddt_df["size"] = plddt_df["seed_id"].str.split("-", expand=True)[2].astype(int) - plddt_df["seed_id"].str.split("-", expand=True)[1].astype(int) + 1
protperties_df = pi_df.merge(ss_info_df, on="seed_id").merge(cn_df, on="seed_id").merge(plddt_df, on="seed_id")




feature_name_mapping = {"avg_contact_num": "Average contact number",
                "avg_intra_fam_pident": "Average sequence identity with other family members",
                "size": "Domain length",
                "avg_plddt" : "Average pLDDT",
                "len_norm_tr_count": "The number of transitions in the secondary structure state normalized by the domain length",
                "c_frac": "Coil fraction",
                "e_frac": "Sheet fraction",
                "h_frac": "Helix fraction"}
tool_name_mapping = {'fs3di': "Foldseek (3Di)", 'fs': "Foldseek (cif_cut)", 'hmmscan': "HMMER", 'mm': "MMseqs2", 'rs':"Reseek", 'tm':"TM-align"}

## Train a model for overall residue alignment
paths = sorted(glob.glob(f"{base_dir}/data/processed/residue_ali_frac_per_seed_split_vs_split/*_all.tsv"))
model_info = {}
for path in paths:
    tool = os.path.basename(path).replace("_all.tsv", "")
    perf_df = pd.read_csv(path, sep="\t").rename(columns={"query": "seed_id"})
    properties_perf = protperties_df.merge(perf_df, on = "seed_id")
    data_type = os.path.basename(path).replace(".tsv", "")
    model_info[data_type] = train_model(properties_perf)
    print(f"trained a model for overall residue alignment of {data_type}")

writer = pd.ExcelWriter(f'{data_dir}/processed/feature_importance_overall_residue_alignment.xlsx', engine='xlsxwriter')

for tool, model_data in model_info.items():
    tool = tool.split("_")[0]
    feature_importance_df = model_data["feature_importance"]
    feature_importance_df["Feature"] = feature_importance_df["Feature"].apply(lambda x: feature_name_mapping[x])
    feature_importance_df.to_excel(writer, sheet_name=tool_name_mapping[tool], index=None)

# Close the Pandas Excel writer and output the Excel file.
writer.close()


## Train a model for conserved residue alignment
paths = sorted(glob.glob(f"{base_dir}/data/processed/residue_ali_frac_per_seed_split_vs_split/*_conserved.tsv"))
model_info = {}
for path in paths:
    tool = os.path.basename(path).replace("_conserved.tsv", "")
    perf_df = pd.read_csv(path, sep="\t").rename(columns={"query": "seed_id"})
    properties_perf = protperties_df.merge(perf_df, on = "seed_id")
    data_type = os.path.basename(path).replace(".tsv", "")
    model_info[data_type] = train_model(properties_perf)
    print(f"trained a model for conserved residue alignment of {data_type}")

writer = pd.ExcelWriter(f'{data_dir}/processed/feature_importance_conserved_residue_alignment.xlsx', engine='xlsxwriter')

for tool, model_data in model_info.items():
    tool = tool.split("_")[0]
    feature_importance_df = model_data["feature_importance"]
    feature_importance_df["Feature"] = feature_importance_df["Feature"].apply(lambda x: feature_name_mapping[x])
    feature_importance_df.to_excel(writer, sheet_name=tool_name_mapping[tool], index=None)

# Close the Pandas Excel writer and output the Excel file.
writer.close()