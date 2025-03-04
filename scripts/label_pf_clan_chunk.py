import pandas as pd
from headers import headers
import os


def label_pf_clan_chunk(input_path,
                        clan_info_path="./data/raw/pfam_clan_info.tsv",
                        remove_cif_ext="Auto",
                        remove_self_matches=True,
                        chunksize=10_000):
    """Requires the input_path. It assumes the input is in the following format: {search_tool}_B{batch_num}.tsv"""
    search_tool = os.path.basename(input_path).split("_B")[0]
    ali_header = headers[search_tool]
    ipr_clan_df = pd.read_csv(clan_info_path, sep="\t")
    if remove_cif_ext == "True":
        remove_cif = True
    elif remove_cif_ext == "False":
        remove_cif = False
    elif remove_cif_ext == "Auto":
        if ".cif" in open(input_path).readline().split()[0]:
            remove_cif = True
        else:
            remove_cif = False

    for chunk in pd.read_csv(input_path, sep="\t", header=None, names=ali_header, chunksize=10_000):
        if remove_cif:
            chunk["query"] =  chunk["query"].str.replace(".cif", "")
            chunk["target"] = chunk["target"].str.replace(".cif", "")
        if remove_self_matches:
            chunk = chunk[chunk["query"] != chunk["target"]]  # This removes self-matches
        chunk["q_pfam"] = chunk["query"].str.split("-", expand=True)[3]
        chunk["t_pfam"] = chunk["target"].str.split("-", expand=True)[3]
        chunk = chunk.reset_index(names=['row_num'])

        chunk = chunk.merge(ipr_clan_df, left_on="q_pfam", right_on="pfam").rename(columns={"clan": "q_clan"}).drop(columns=["pfam"])
        chunk = chunk.merge(ipr_clan_df, left_on="t_pfam", right_on="pfam").rename(columns={"clan": "t_clan"}).drop(columns=["pfam"])
        chunk["pfam_label"] = (chunk["q_pfam"] == chunk["t_pfam"])
        chunk["clan_label"] = (chunk["q_clan"] == chunk["t_clan"])

        yield chunk
