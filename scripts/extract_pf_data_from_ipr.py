from pathlib import Path
import os
import pandas as pd


try:
    # script location -> parent of the script's directory (project root)
    base_dir = Path(__file__).resolve().parent.parent
except NameError:
    # running in a notebook or interactive shell where __file__ is not defined
    base_dir = Path.cwd().resolve().parent


dom_data_df = pd.read_csv(f"{base_dir}/data/tb_ipr.tsv", sep="\t", header=None)
pf_data = dom_data_df[dom_data_df[3].str.startswith("PF")][[0, 3, 4, 5]].reset_index(drop=True)
pf_data.to_csv(f"{base_dir}/data/tb_pfam.tsv", sep="\t", header=None, index=None)
