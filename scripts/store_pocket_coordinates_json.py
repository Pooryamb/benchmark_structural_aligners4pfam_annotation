import glob
from concurrent.futures import ProcessPoolExecutor
import json
import warnings
import pandas as pd

warnings.filterwarnings("ignore")

prank_preds = glob.glob("./tmp/prank_preds/*.cif_predictions.csv")
conf_threshold = 0.8
path_df = pd.DataFrame({"path": prank_preds})

path_df["family"] = path_df["path"].str.split("-", expand=True)[3].str.split(".", expand=True)[0]
path_df["seed_id"] = path_df["path"].str.split("/", expand=True)[3].str.replace(".cif_predictions.csv", "")

subdfs = [group for _, group in path_df.groupby("family")] # By this end, we have a dataframe for each family containing its information

def process_fam_members(fam_df):
    pockets_residues = {}
    family_id = fam_df.iloc[0]["family"]
    for index, row in fam_df.iterrows():
        seed_df = pd.read_csv(row["path"], delimiter=",", skipinitialspace=True)
        high_conf_pockets = seed_df[seed_df["probability"]>=conf_threshold]
        high_conf_pockets["residue_ids"] = high_conf_pockets["residue_ids"].str.replace("A_", "")
        all_high_conf_pocket_residues = [int(x) for x in " ".join(high_conf_pockets["residue_ids"]).split()]
        pockets_residues[row["seed_id"]] = all_high_conf_pocket_residues
    return pockets_residues

with ProcessPoolExecutor() as executor:
    results = executor.map(process_fam_members, subdfs)

single_dict_result = {}
for member_dict in results:
    single_dict_result.update(member_dict)

with open("./tmp/residue_features/pocket.json", 'w') as ofile:
    json.dump(single_dict_result, ofile, indent=4)
