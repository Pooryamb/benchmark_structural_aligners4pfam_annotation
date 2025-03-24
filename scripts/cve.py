import os
import pandas as pd
import numpy as np


def cve(df):
    """Takes a dataframe whose rows show for each query within each e-value evalue_bin
    how many True Positives or False Positives are found. It returns a dataframe
    that for each e-value threshold has the coverage and fpepq.
    The input dataframe should have the following columns:
        query: The accession of the query
        evalue_bin: This is the round(log of e-value)
        tp_pfam: The number of TPs at Pfam level at this bin
        tp_clan: The number of TPs at Clan level at this bin
        fp_pfam: The number of FPs at Pfam level at this bin
        fp_clan: The number of FPs at Clan level at this bin
    """
    cnt_cols = ["tp_pfam", "tp_clan", "fp_pfam", "fp_clan"]
    res_df = df.groupby("evalue_bin")[cnt_cols].agg(sum).reset_index().sort_values(by="evalue_bin")
    q_count = df["query"].unique().shape[0]
    #As FEPEQ is False positives expected per query, we need to divide the number of FPs by the number of queries.
    res_df[["fp_pfam", "fp_clan"]] = res_df[["fp_pfam", "fp_clan"]]/q_count
    res_df["fpepq_pfam"] = res_df["fp_pfam"].cumsum()
    res_df["fpepq_clan"] = res_df["fp_clan"].cumsum()
    res_df["cov_pfam"] = res_df["tp_pfam"].cumsum()/res_df["tp_pfam"].sum()
    res_df["cov_clan"] = res_df["tp_clan"].cumsum()/res_df["tp_clan"].sum()
    return res_df
