import glob
import os
import pandas as pd
import numpy as np
import xlsxwriter
from sklearn.model_selection import train_test_split  
from sklearn.ensemble import RandomForestRegressor  
from sklearn.metrics import mean_squared_error
from sens_up2_first_fp import get_cum_sens_from_1st_occ_df, get_tps_frac_bef_1st_fp_from_1st_occ_df, cum_sens_dist, auc_tps_bef_1st_fp

script_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(script_dir, "../")
data_dir = f"{base_dir}/data/"

# Load domain features
pi_df = pd.read_csv(f"{data_dir}/processed/avg_intra_fam_pident.tsv", sep="\t") #pi means percentage identity
ss_info_df = pd.read_csv(f"{data_dir}/processed/ss_info_pfam.tsv", sep="\t")
cn_df = pd.read_csv(f"{data_dir}/processed/avg_contact_num.tsv", sep="\t")
plddt_df = pd.read_csv(f"{data_dir}/processed/pfam_avg_plddt.tsv", sep="\t")
plddt_df["size"] = plddt_df["seed_id"].str.split("-", expand=True)[2].astype(int) - plddt_df["seed_id"].str.split("-", expand=True)[1].astype(int) + 1

protperties_df = pi_df.merge(ss_info_df, on="seed_id").merge(cn_df, on="seed_id").merge(plddt_df, on="seed_id")

# Find SFFP for each query

first_occ_dict = {}

all_tsv_files = glob.glob(f"{data_dir}/processed/first_label_occ/*_B*.tsv")
tools = {'cif_cut', 'mm', 'reseek', 'tm'} 
file_paths_tool_dict = {x: [y for y in all_tsv_files if x in y] for x in tools} 

for tool, tool_tsv_paths in file_paths_tool_dict.items():
    first_occ_list = []
    for tsv_path in tool_tsv_paths:
        tsv_df = pd.read_csv(tsv_path, sep="\t")
        first_occ_list.append(tsv_df)
    first_occ_dict[tool] = pd.concat(first_occ_list)
    

frac_sens_up2first_fp_dict_fam = {tool: get_tps_frac_bef_1st_fp_from_1st_occ_df(first_occ_dict[tool])[["seed_id", "tp_bef_fp_frac_pfam"]] for tool in first_occ_dict.keys()}
perf_prop = {tool: frac_sens_up2first_fp_dict_fam[tool].merge(protperties_df, on="seed_id") for tool in tools}

# A function to train a model to predict performance based on other features

def train_model(data_df):
    """Trains a random forest for each model and shows the importance of the features"""
    # Assuming the last 8 columns are features and the second column is the target  
    X = data_df.iloc[:, -8:]     # Features (last 8 columns)  
    y = data_df.iloc[:, 1]       # Target (second column)
    
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
    

tool_information = {}
for tool, data_df in perf_prop.items():
    tool_information[tool] = train_model(data_df)

feature_name_mapping = {"avg_contact_num": "Average contact number",
                "avg_intra_fam_pident": "Average sequence identity with other family members",
                "size": "Domain length",
                "avg_plddt" : "Average pLDDT",
                "len_norm_tr_count": "The number of transitions in the secondary structure state normalized by the domain length",
                "c_frac": "Coil fraction",
                "e_frac": "Sheet fraction",
                "h_frac": "Helix fraction"}
tool_name_mapping = {'cif_cut': "Foldseek (cif_cut)", 'mm': "MMseqs", 'reseek': "Reseek", 'tm': "TM-align"} 


# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter(f'{data_dir}/processed/feature_importance_sffp.xlsx', engine='xlsxwriter')


for tool, model_data in tool_information.items():
    feature_importance_df = model_data["feature_importance"]
    feature_importance_df["Feature"] = feature_importance_df["Feature"].apply(lambda x: feature_name_mapping[x])
    feature_importance_df.to_excel(writer, sheet_name=tool_name_mapping[tool], index=None)

# Close the Pandas Excel writer and output the Excel file.
writer.close()