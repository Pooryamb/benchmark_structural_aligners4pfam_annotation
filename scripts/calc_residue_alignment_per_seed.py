import glob
import os
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

def summarize_family(tsv_path):
    res_ali_df = pd.read_csv(tsv_path, sep="\t")
    res_ali_df = res_ali_df[res_ali_df["all"]!=0]
    res_ali_df = res_ali_df[res_ali_df["query"]!=res_ali_df["target"]]
    res_ali_df["aligned_fraction"] = res_ali_df["correct"] / res_ali_df["all"]
    return res_ali_df.groupby("target")["aligned_fraction"].agg("mean").reset_index()


tools = [os.path.basename(x) for x in glob.glob("./tmp/intrafam_residue_alignment_counts/*")]
res_types = [os.path.basename(x) for x in glob.glob(f"./tmp/intrafam_residue_alignment_counts/{tools[0]}/*")]


for ali_tool in tools:
    for kind in res_types:
        tsv_paths = glob.glob(f"./tmp/intrafam_residue_alignment_counts/{ali_tool}/{kind}/*.tsv")
        with ProcessPoolExecutor() as executor:
            results = executor.map(summarize_family, tsv_paths)
            results = [df for df in results if df.shape[0]>0]
            try:
                results_df = pd.concat(results)
                results_df.to_csv(f"./data/processed/residue_ali_frac_per_seed/{ali_tool}_{kind}.tsv", sep="\t", index=None)
            except:
                print(ali_tool)
                print(kind)
