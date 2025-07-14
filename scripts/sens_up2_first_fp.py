import glob
import os
import pandas as pd
import numpy as np


script_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(script_dir, "../data")

def cnt_tps_bef_1st_fp(input_df):
    """Takes a dataframe containing the first occurence of each (pfam_label and clan_label) as input
    and returns a dataframe, where each row has a query name, the number of tps before the
    first fp at pfam level and clan level.
    """
    tps_bef_first_fp = []
    grouped = input_df.groupby("query")
    for seed_id, seed_info in grouped:
        min_rnum = seed_info.head(1)["row_num"]
        first_fp_pfam_rnum = seed_info[seed_info["pfam_label"]==False].drop_duplicates(["pfam_label"])["row_num"]
        first_fp_clan_rnum = seed_info[seed_info["clan_label"]==False].drop_duplicates(["clan_label"])["row_num"]
        tps_bef_first_fp_pfam = list(first_fp_pfam_rnum)[0] - list(min_rnum)[0]
        tps_bef_first_fp_clan = list(first_fp_clan_rnum)[0] - list(min_rnum)[0]
        tps_bef_first_fp.append([seed_id, tps_bef_first_fp_pfam, tps_bef_first_fp_clan])
    return pd.DataFrame(tps_bef_first_fp, columns=["seed_id", "pfam_tp", "clan_tp"])


#Now, we find the number of seeds in each Pfam/Clan
ipr_data = pd.read_csv(f"{data_dir}/raw/pfam_clan_info.tsv", sep="\t")
fam_size_df = pd.read_csv(f"{data_dir}/raw/dbs/pfam_fs_cut_clust/pfam_h.tsv", sep="\t", header=None)
fam_size_df[ "pfam"] = fam_size_df[1].str.split("-", expand=True)[3]
fam_size_df = fam_size_df.merge(ipr_data, on="pfam").drop(columns=[0])
fam_size_df["cnt_pfam"] = fam_size_df.groupby("pfam")["pfam"].transform("count")
fam_size_df["cnt_clan"] = fam_size_df.groupby("clan")["clan"].transform("count")
fam_size_df = fam_size_df.drop(columns=["pfam", "clan"])
fam_size_df = fam_size_df.rename(columns={1: "seed_id"})


def calc_tps_frac_bef_1st_fp(cnt_df, fam_size_df, ignore_self_hits=True):
    """Requires two dataframes as inpu:
       dataframe 1: the number of tps before first fp at pfam and clan level
       dataframe 2: the number of seeds in each family/clan
       it is possible to specify if self_hits has to be ignored. By default, they will be ignored.
       The first dataframe can be provided by cnt_tps_bef_1st_fp function.
       The second dataframe can be calculated by off-function scripts.
       """
    freq_df = cnt_df.merge(fam_size_df, on="seed_id")
    deduct = 0
    if ignore_self_hits:
        deduct = 1
        freq_df = freq_df[freq_df["cnt_pfam"]>1]
    freq_df["tp_bef_fp_frac_pfam"] = freq_df["pfam_tp"] / (freq_df["cnt_pfam"] - deduct)
    freq_df["tp_bef_fp_frac_clan"] = freq_df["clan_tp"] / (freq_df["cnt_clan"] - deduct)
    freq_df = freq_df.reset_index(drop=True)
    return freq_df


def cum_sens_dist(tps_frac, target_col="tp_bef_fp_frac_pfam"):
    """ The inputs are a dataframe containing the fraction of tps before first fp and the column to calculate
    cumulative frequency based on. The input dataframe can be provided by calc_tps_frac_bef_1st_fp function
    """
    sorted_df = tps_frac.sort_values(by = target_col, ascending=False)
    sorted_df["row_num"] = 1
    sorted_df["cum_sum"] = sorted_df["row_num"].cumsum()
    sorted_df["cum_sum_frac"] = sorted_df["cum_sum"]/tps_frac.shape[0]
    sorted_df = sorted_df.drop(columns=["row_num", "cum_sum"])
    return sorted_df

def get_tps_frac_bef_1st_fp_from_1st_occ_df(first_occ_df):
    cnt_tps = cnt_tps_bef_1st_fp(first_occ_df)
    frac_tps = calc_tps_frac_bef_1st_fp(cnt_tps, fam_size_df)
    return frac_tps

def get_cum_sens_from_1st_occ_df(first_occ_df, target_col):
    """Takes the dataframe as input and returns the fraction of seeds with a sensivity up to a threshold as output"""
#    cnt_tps = cnt_tps_bef_1st_fp(first_occ_df)
    frac_tps = get_tps_frac_bef_1st_fp_from_1st_occ_df(first_occ_df)
    #calc_tps_frac_bef_1st_fp(cnt_tps, fam_size_df)
    return cum_sens_dist(frac_tps, target_col)

def auc_tps_bef_1st_fp(input_df, target_col="tp_bef_fp_frac_pfam"):
    """This function takes a dataframe containing the fraction of tps before first fp as input
    and returns the area under the cumulative distribution of fraction of tps before first fp as the output.
    """
    cum_df = cum_sens_dist(input_df, target_col)
    x_axis = cum_df["cum_sum_frac"]
    y_axis = cum_df[target_col]
    auc = np.trapz(y_axis, x_axis)
    return auc
