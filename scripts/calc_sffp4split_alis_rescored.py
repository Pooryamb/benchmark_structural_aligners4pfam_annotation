import glob
import sys
import os
import multiprocessing as mp
import pickle
import pandas as pd


script_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(script_dir, "..")

ipr_data = pd.read_csv(f"{base_dir}/data/raw/pfam_clan_info.tsv", sep="\t")

# Find the number of members in pfam, clan
fam_size_df = pd.read_csv(f"{base_dir}/data/raw/dbs/pfam_split_target/pfam_h.tsv", sep="\t", header=None)
fam_size_df["pfam"] = fam_size_df[1].str.split("-", expand=True)[3]
fam_size_df = fam_size_df.merge(ipr_data, on="pfam").drop(columns=[0])
fam_size_df["pfam_size"] = fam_size_df.groupby("pfam")["pfam"].transform("count")
fam_size_df["clan_size"] = fam_size_df.groupby("clan")["clan"].transform("count")
fam_size_df = fam_size_df.drop_duplicates(["pfam", "clan", "pfam_size", "clan_size"]).drop(columns=[1, "clan"]).reset_index(drop=True)

def read_and_prepare_ali_data(input_path, cols2read):
    """Reads alignment data and labels each row"""
    ali_df = pd.read_csv(input_path, sep="\t", usecols=cols2read)
    # Add q_pfam data
    ali_df["q_pfam"] = ali_df["query"].str.split("-", expand=True)[3]
    ali_df = ali_df.merge(ipr_data, left_on="q_pfam", right_on="pfam").drop(columns=["pfam"]).rename(columns={"clan": "q_clan"})
    # Add t_pfam data
    ali_df["t_pfam"] = ali_df["target"].str.split("-", expand=True)[3]
    ali_df = ali_df.merge(ipr_data, left_on="t_pfam", right_on="pfam").drop(columns=["pfam"]).rename(columns={"clan": "t_clan"})
    # Add pfam label
    ali_df["pfam_label"] = (ali_df["t_pfam"] == ali_df["q_pfam"])
    # Add clan label
    ali_df["clan_label"] = (ali_df["t_clan"] == ali_df["q_clan"])
    return ali_df

def select_first_occ_each_label(ali_df, col2sortby = "evalue", sort_ascending=True):
    """Takes the read and labeled dataset as input, sorts based on the column specified, and then selects the first occurence of each label at Pfam/clan level"""
    ali_df = ali_df.sort_values(by=["query", col2sortby], ascending=[True, sort_ascending], ignore_index=True) # evalue can be replaced by other things such as alignment score (which also changes the ascending to descending)
    ali_df = ali_df.reset_index(names=["row_num"])
    first_occ_df = ali_df.drop_duplicates(["query", "pfam_label", "clan_label"]) # Stores the first occurence of each label set
    return first_occ_df


def count_tps_before_first_fp_per_query(grp_df):
    """Receving a dataframe containing the row_number for the first TP/FP at clan/pfam level, it returns the number of TPs before the first
    FP at Pfam/Clan level"""
    min_row_num = grp_df.iloc[0]["row_num"] # The index of the first row

    fps_pfam = grp_df[~(grp_df["pfam_label"])]
    if len(fps_pfam) == 0:
        first_fp_pfam_rownum = min_row_num
    else:
        first_fp_pfam_rownum = fps_pfam.iloc[0]["row_num"]

    fps_clan = grp_df[~(grp_df["clan_label"])]
    if len(fps_clan) == 0:
        first_fp_clan_rownum = min_row_num
    else:
        first_fp_clan_rownum = fps_clan.iloc[0]["row_num"]

    tps_bef_1st_fp_pfam = first_fp_pfam_rownum - min_row_num  # If the first FP shows before the first TP, it must return 0
    tps_bef_1st_fp_clan = first_fp_clan_rownum - min_row_num  # If the first FP shows before the first TP, it must return 0

    return tps_bef_1st_fp_pfam, tps_bef_1st_fp_clan

def count_tps_before_first_for_first_occ_df(first_occ_df):
    """Receiving a dataframe of the first occurences, it will return a dataframe containing the number of TPs/FPs for each query at Pfam/Clan level"""
    count_up2first_fp = first_occ_df.groupby("query")[["row_num", "pfam_label", "clan_label"]].apply(count_tps_before_first_fp_per_query)
    count_up2first_fp = pd.DataFrame(count_up2first_fp.tolist(), index=count_up2first_fp.index, columns=['tps_bef_first_fp_pfam', 'tps_bef_first_fp_clan']).reset_index()
    return count_up2first_fp


def process_one_file_path_to_find_tps_count_bef_first_fp(input_path, col2sortby, sort_ascending, cols2read):
    """Encapsulated the operations for a single file."""
    print(f"Processing {input_path}")
    ali_df = read_and_prepare_ali_data(input_path, cols2read)

    # Calculate TPs and FPs
    max_tps_pfam = ali_df.groupby("query")["pfam_label"].sum().reset_index().rename(columns={"pfam_label": "max_tps_pfam"})
    max_tps_clan = ali_df.groupby("query")["clan_label"].sum().reset_index().rename(columns={"clan_label": "max_tps_clan"})

    first_occ_df = select_first_occ_each_label(ali_df, col2sortby, sort_ascending)

    count_up2first_fp = count_tps_before_first_for_first_occ_df(first_occ_df)

    # Merging results
    tps_bef_1st_fp_data = count_up2first_fp.merge(max_tps_pfam, on="query").merge(max_tps_clan, on="query")
    tps_bef_1st_fp_data["pfam"] = tps_bef_1st_fp_data["query"].str.split("-", expand=True)[3]

    # Assuming fam_size_df is globally accessible or you can pass it as a parameter
    tps_bef_1st_fp_data = tps_bef_1st_fp_data.merge(fam_size_df, on="pfam").drop(columns=["pfam"])
    return tps_bef_1st_fp_data


def find_tps_count_bef_first_fp_for_filepaths(filepaths, cols2read, col2sortby="evalue", sort_ascending=True):
    """Receive the filepaths and calculate the number of TPs before the first FP for all files and return a single DataFrame as output."""

    # Parallel processing of file paths ##Uses one fourth of the data because of RAM demanding operation
    with mp.Pool(4) as pool: #processes=mp.cpu_count()
        # Use lambda to inject additional parameters into the process_one_file_path_to_find_tps_count_bef_first_fp function
        batch_data_list = pool.starmap(process_one_file_path_to_find_tps_count_bef_first_fp, [(path, col2sortby, sort_ascending, cols2read) for path in filepaths])

    # Concatenating the results
    concat_df = pd.concat(batch_data_list)
    concat_df["sffp_pfam"] = concat_df["tps_bef_first_fp_pfam"] / concat_df["pfam_size"]
    concat_df["sffp_clan"] = concat_df["tps_bef_first_fp_clan"] / concat_df["clan_size"]
    concat_df["max_sffp_pfam"] = concat_df["max_tps_pfam"] / concat_df["pfam_size"]
    concat_df["max_sffp_clan"] = concat_df["max_tps_clan"] / concat_df["clan_size"]

    concat_df = concat_df.reset_index(drop=True)
    return concat_df



foldseek_paths = glob.glob(f"{base_dir}/tmp/alis/split_pf_seq/fs_pref*_bitscore.tsv")
reseek_paths = glob.glob(f"{base_dir}/tmp/alis/split_pf_seq/reseek_sens_*_bitscore.tsv")
cols2read = ["query", "target", "evalue", "pssm_raw_score", "bitscore_rep"]
tps_bef_first_fp = {
                    "foldseek_bitscore": find_tps_count_bef_first_fp_for_filepaths(foldseek_paths, cols2read, col2sortby="bitscore_rep", sort_ascending=False),
                    "foldseek": find_tps_count_bef_first_fp_for_filepaths(foldseek_paths,cols2read, col2sortby="evalue", sort_ascending=True),
                    "foldseek_rescored": find_tps_count_bef_first_fp_for_filepaths(foldseek_paths,cols2read, col2sortby="pssm_raw_score", sort_ascending=False),
                    "reseek": find_tps_count_bef_first_fp_for_filepaths(reseek_paths, cols2read,col2sortby="evalue", sort_ascending=True),
                    "reseek_rescored": find_tps_count_bef_first_fp_for_filepaths(reseek_paths, cols2read, col2sortby="pssm_raw_score", sort_ascending=False),
                    "reseek_bitscore": find_tps_count_bef_first_fp_for_filepaths(reseek_paths, cols2read, col2sortby="bitscore_rep", sort_ascending=False)
                    }

sffp_dir_path = f"{base_dir}/tmp/sffp_rescoring"
os.makedirs(sffp_dir_path, exist_ok=True)
for key, value in tps_bef_first_fp.items():
    value.to_csv(f"{sffp_dir_path}/{key}.tsv", sep="\t", index=None)
