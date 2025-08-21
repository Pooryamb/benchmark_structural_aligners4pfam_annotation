import pickle
import os
import sys
import glob
import pandas as pd
import numpy as np
import numba as nb
import swifter
from headers import fam_ali_headers
from concurrent.futures import ProcessPoolExecutor
import argparse


script_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(script_dir, "..")


parser = argparse.ArgumentParser(description="Calculates the psi-blast raw score for the alignments provided either by Foldseek or Reseek")
parser.add_argument("--alis", type=str, required=True, help="This could be r (for all reseek alignments), f (for all foldseek alignments, a for all foldseek and reseek, or a filepath for a specific file")
args = parser.parse_args()

if args.alis == "f":
    inp_paths = glob.glob(f"{base_dir}/tmp/alis/split_pf_seq/fs_pref*")
elif args.alis == "r":
    inp_paths = glob.glob(f"{base_dir}/tmp/alis/split_pf_seq/reseek_sens*")
elif args.alis == "a":
    inp_paths = glob.glob(f"{base_dir}/tmp/alis/split_pf_seq/fs_pref*") + glob.glob(f"{base_dir}/tmp/alis/split_pf_seq/reseek_sens*")
else:
    inp_paths = [args.alis]

inp_paths = [x for x in inp_paths if "rescored" not in x]

@nb.njit
def map_residues_with_gaps_numba(qseq, tseq, qstart=1, tstart=1):
    n = len(qseq)
    mapping = np.zeros((n, 2), dtype=np.int16)
    qpos, tpos = qstart, tstart
    for i in range(n):
        qc = qseq[i]
        tc = tseq[i]
        if qc != ord('-') and tc != ord('-'):
            mapping[i, 0] = qpos
            mapping[i, 1] = tpos
            qpos += 1
            tpos += 1
        elif qc != ord('-'):
            mapping[i, 0] = qpos
            mapping[i, 1] = -1
            qpos += 1
        elif tc != ord('-'):
            mapping[i, 0] = -1
            mapping[i, 1] = tpos
            tpos += 1
    return mapping


def map_row_residues(row):
    # Convert string to bytes (Numba works with bytes, not Python str) It is for speed up
    qaln = row.qaln.encode()
    taln = row.taln.encode()
    return pd.DataFrame(map_residues_with_gaps_numba(qaln, taln, row.qstart, row.tstart), columns=["q_coord", "seed_coord"])

def find_gap_openings_and_extensions(gap_df, column_wo_gap):
    """Finds the number of gap opening and gap extensions, assuming that the rows are sorted. It also requires the column which
    shows the entry without gaps. """
    gap_openings = (gap_df[column_wo_gap].diff()!=1).sum()
    gap_extensions = len(gap_df) - gap_openings
    return gap_openings, gap_extensions

def calc_gap_score(gap_info):
    """given the gap information, it calculates the gap score based on the number of gap openings and gap extensions"""
    gap_opening_penalty = -11
    gap_extension_penalty = -1
    gap_penalty = gap_info[0] * gap_opening_penalty + gap_info[1] * gap_extension_penalty
    return gap_penalty

def get_ali_score_for_row(row, pssms, seed2cons_mapping):
    """rescores a row based on pssm and residue mapping of the target seed to the consensus sequence"""
    mapping = map_row_residues(row) # Step 1: get mapping between q residues and seed residues
    mapping["q_residue"] = list(row["qaln"])       # Add q sequence
    q_res_is_gap = mapping["q_coord"] == -1        # Index of rows that q has a gap
    seed_res_is_gap = mapping["seed_coord"] == -1  # Index of rows that seed has a gap
    non_gap_part = mapping[~( q_res_is_gap| seed_res_is_gap)]
    non_gap_part = non_gap_part.merge(seed2cons_mapping[row["target"]], on="seed_coord")
    non_gap_part = non_gap_part.merge(aa2ind, left_on="q_residue", right_on="aa")
    pssm = pssms[row["t_family"]]                  # Access the PSSM of the family
    pssm_rows2sel = (non_gap_part["cons_coord"] - 1).to_numpy() # Proper rows of PSSM (PSSM row numbers start from 0, whereas mappings start from 1)
    pssm_cols2sel = non_gap_part["aa_ind"].to_numpy()           # Proper columns of PSSM
    no_gap_scores = pssm.values[pssm_rows2sel, pssm_cols2sel].sum()      # Select proper values from PSSM

    gap_on_q = mapping[q_res_is_gap]               # Rows that q has a gap
    q_gap_penalty = calc_gap_score(find_gap_openings_and_extensions(gap_on_q, "seed_coord")) # The total alignment score for the ungapped part

    gap_on_seed = mapping[seed_res_is_gap]         # Rows that seed has a gap
    seed_gap_penalty = calc_gap_score(find_gap_openings_and_extensions(gap_on_seed, "q_coord")) # The total alignment score for the ungapped part

    return no_gap_scores + q_gap_penalty + seed_gap_penalty


def split_scoring_data(ali_df, pssms, seed2cons_mapping):

    """It splits the alignment dataframe, PSSMs, and residue mappings into multiple chunks for speed up"""
    ali_df = ali_df.sort_values(by="t_family").reset_index(drop=True)
    cpu_count = os.cpu_count()
    ali_dfs = np.array_split(ali_df, cpu_count)
    chunked_data = []
    for chunk in ali_dfs:
        chunk_t_family_list = chunk["t_family"].unique()
        chunk_targets = chunk["target"].unique()
        chunk_pssms = {x:pssms[x] for x in chunk_t_family_list}
        chunk_seed2cons_mapping = {x: seed2cons_mapping[x] for x in chunk_targets}
        chunked_data.append({"ali_df": chunk, "pssms": chunk_pssms, "seed2cons_mapping": chunk_seed2cons_mapping})
    return chunked_data


def rescore_dataframe(ali_df, pssms, seed2cons_mapping):
    """It will rescore a dataframe based on provided PSSMs and residue mappings"""
    ali_df["pssm_raw_score"] = ali_df.apply(lambda row: get_ali_score_for_row(row, pssms, seed2cons_mapping), axis=1)  # .apply(get_ali_score, axis=1)
    ali_df = ali_df.drop(["qaln", "taln"], axis=1)
    return ali_df


def rescore_handler(args):
    """Defined to facilitate running the scoring function in parallel mode"""
    return rescore_dataframe(args["ali_df"], args["pssms"], args["seed2cons_mapping"])


# Load the PSSM data
with open(f"{base_dir}/data/processed/parsed_pssm.pkl", 'rb') as file:
    pssm_data_dict = pickle.load(file)

pssms = pssm_data_dict["pssm_profiles"]
aa2ind = pd.DataFrame(pssm_data_dict["aa_ind"].items(), columns=["aa", "aa_ind"])
k_values =  pssm_data_dict["k_values"]
lamda_values =  pssm_data_dict["lamda_value"]
# Load residue mapping between target and the consensus sequence
with open(f"{base_dir}/data/processed/seed2cons_mapping.pkl", 'rb') as file:
    seed2cons_mapping = pickle.load(file)

fs_header = fam_ali_headers["fs"][:]
fs_header.remove("lddtfull")

rs_header = fam_ali_headers["rs"]


#inp_paths = glob.glob(f"{base_dir}/tmp/alis/split_pf_seq/fs_pref*") + glob.glob(f"{base_dir}/tmp/alis/split_pf_seq/reseek_sens*")
#inp_paths = [x for x in inp_paths if "rescored" not in x]

for inp_path in inp_paths:
    print(f"processing {inp_path}")
    if "reseek" in os.path.basename(inp_path):
        columns = rs_header
    elif "fs_pref" in os.path.basename(inp_path):
        columns = fs_header
    ali_df = pd.read_csv(inp_path, sep="\t", header=None, names=columns, keep_default_na=False, na_values=[""])
    ali_df["t_family"] = ali_df["target"].str.split("-", expand=True)[3]
    split_chunks = split_scoring_data(ali_df, pssms, seed2cons_mapping)
# Using a single CPU:
#    rescored_chunks = []
#    for chunk in split_chunks:
#        rescored_chunk = rescore_handler(chunk)
#        rescored_chunks.append(rescored_chunk)

    with ProcessPoolExecutor() as executor:
    # Map rescore_handler over split_chunks in parallel
        rescored_chunks = list(executor.map(rescore_handler, split_chunks))
    all_alidf_rescored = pd.concat(rescored_chunks)
    all_alidf_rescored.to_csv(inp_path.replace(".tsv", "_rescored.tsv"), sep="\t", index=None)
