import os
import pandas as pd


script_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(script_dir, "..")

df = pd.read_csv(f"{base_dir}/data/raw/dbs/tb/tb.lookup", sep="\t", header=None)
unip_ids = df[1].str.split("-", expand=True)[1]
unip_ids = unip_ids.sort_values().reset_index(drop=True)
unip_ids.to_csv(f"{base_dir}/data/tb_unip_ids.tsv", sep="\t", header=None, index=None)

