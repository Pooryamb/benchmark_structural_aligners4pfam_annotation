import glob
import multiprocessing as mp
import os
import pandas as pd


script_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(script_dir, "..")

#base_dir = "."

rescored_files = glob.glob(f"{base_dir}/tmp/alis/split_pf_seq/*_rescored.tsv")
out_dir = f"{base_dir}/tmp/structurally_similar_pairs/"

os.makedirs(out_dir, exist_ok=True)
cols2read = ["query", "target", "evalue"]

def select_equal_num_of_tps_fps(file_path):
    """Sorts hits by evalue and finds the largest evalue where hits up to that evalue are equally TP or FP"""
    df = pd.read_csv(file_path, sep="\t", usecols=cols2read)
    df["qlen"] = df["query"].str.split("-", expand=True)[2].astype(int) - df["query"].str.split("-", expand=True)[1].astype(int) + 1
    df = df[df["qlen"] >= 50] # Only select the queries that are long enough
    df = df.sort_values(by="evalue", ignore_index=True)
    df["label"] = (df["query"].str.split("-", expand=True)[3] == df["target"].str.split("-", expand=True)[3])
    df["cum_tp"] = df["label"].cumsum()
    df["cum_fp"] = df.index - df["cum_tp"] + 1
    df["more_tp"] = (df["cum_tp"] >= df["cum_fp"])
    last_index = df[df['more_tp']].last_valid_index()
    selected_subdf = df.loc[:last_index]
    file_name = os.path.basename(file_path)
    selected_subdf[["query", "target"]].to_csv(f"{out_dir}/{file_name}", sep="\t", header=None, index=None)


with mp.Pool() as pool: #processes=mp.cpu_count()
    pool.map(select_equal_num_of_tps_fps, rescored_files)
