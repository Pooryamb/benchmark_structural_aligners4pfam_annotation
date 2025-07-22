import os
import glob
import json
import pandas as pd


script_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(script_dir, "../data")
base_dir = os.path.join(script_dir, "../")


all_plddts = json.loads(open(f"{data_dir}/pfam_plddt.json").read())
conserved_residues = json.loads(open(f"{base_dir}/tmp/residue_features/conserved.json").read())

def avg(num_list):
    return sum(num_list)/len(num_list)

cons_minus_bg_plddt = {}
for seed_id, con_res in conserved_residues.items():
    if len(con_res) == 0:
        continue
    bkg_plddt = all_plddts[seed_id]
    bkg_avg_plddt = avg(bkg_plddt)
    conserved_residues_avg_plddt = avg([bkg_plddt[i-1] for i in con_res])
    cons_minus_bg_plddt[seed_id] = conserved_residues_avg_plddt - bkg_avg_plddt

data = pd.DataFrame({ "seed_id": cons_minus_bg_plddt.keys(), "conserved_residue_plddt_minus_bkg": cons_minus_bg_plddt.values()})
sorted_data = data.sort_values(by="conserved_residue_plddt_minus_bkg" ,ascending = False).reset_index(drop=True)

more_confident_than_background = (sorted_data["conserved_residue_plddt_minus_bkg"] > 0).argmin() / len(sorted_data)
print(f"For {round(100 * more_confident_than_background, 2)} percent of the seeds, the pLDDT of the conserved sites is higher than the rest of the seed")

seeds_with_low_confidence_conserved_sites = sorted_data[sorted_data["conserved_residue_plddt_minus_bkg"]<0]["seed_id"]

seeds_with_low_confidence_conserved_sites.to_csv(f"{data_dir}/processed/low_confidence_conserved_residues.tsv", sep="\t", index=None)