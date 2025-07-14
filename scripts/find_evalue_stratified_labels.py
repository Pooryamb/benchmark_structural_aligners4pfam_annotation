import argparse
import os
import pandas as pd
import numpy as np
from headers import headers
from label_pf_clan_chunk import label_pf_clan_chunk

parser = argparse.ArgumentParser(description="This script takes "
    "the output alignment of a tool as input and outputs the number of TPs "
    "found at each e-value threshold. It also requires "
    "a file showing the associations between Pfam ID and Clan ID.")

parser.add_argument("--input", type=str, required=True, help="The path to the input alignment file")
parser.add_argument("--output", type=str, help="The number of TPs found for each query at each "
"evalue threshold")
parser.add_argument("--pfam_clan_info", default="./data/raw/pfam_clan_info.tsv", type=str, help= "The path to the file relating pfam to clan")
parser.add_argument("--remove_cif_ext", type=str, default="Auto", help= "Do you want the cif extension? True/False/Auto")

args = parser.parse_args()

if not args.output:
    args.output = "data/processed/evalue_bins/" + os.path.basename(args.input)
os.system(f"rm -f {args.output}")

def find_tps_freq_in_batch(chunk):
    grouped = chunk.groupby(['query', 'evalue_bin'])
    # Calculate True Positives (sum of True values)
    tp = grouped[['pfam_label', 'clan_label']].sum()
    # Calculate Total Occurrences (count of True + False)
    total = grouped[['pfam_label', 'clan_label']].count()
    # Calculate False Positives (total - true positives)
    fp = total - tp
    tp.columns = ['tp_pfam', 'tp_clan']
    fp.columns = ['fp_pfam', 'fp_clan']
    # Combine True Positives and False Positives into a single DataFrame
    result = pd.concat([tp, fp], axis=1).reset_index()
    return result

buffer_df = pd.DataFrame()
processed_df = pd.DataFrame()
last_id = ""
q_id = 0
is_first_q = False
is_tm = "tm_" in args.input

for chunk in label_pf_clan_chunk(args.input, remove_self_matches=False):
#Please note that the chunks should not go above 250k which is the size of the target database. It is 10k by default
    chunk["evalue"] = pd.to_numeric(chunk["evalue"], errors="coerce")
    if is_tm:
        chunk["evalue_bin"] = -np.round(100 * chunk['evalue'])
    else:
        chunk["evalue_bin"] = np.where(chunk['evalue'] == 0, -1000, np.round(np.log2(chunk['evalue'])))
    head_id = chunk.iloc[0]["query"]
    if not(buffer_df.empty):
        chunk = pd.concat([buffer_df, chunk])
        buffer_df = pd.DataFrame()
    if head_id != last_id and last_id != "":
        if q_id == 0:
            is_first_q = True
            q_id += 1
        else:
            is_first_q = False
        processed_df.reset_index().to_csv(args.output, sep="\t", header=is_first_q, index=False, mode='a')
        processed_df = pd.DataFrame()

    head_chunk = chunk[chunk["query"] == head_id]
    # Handle addition of new batch to the processed DataFrame
    head_chunk_result = find_tps_freq_in_batch(head_chunk).set_index(['query', 'evalue_bin'])
    if processed_df.empty:
        # If processed_df is empty, directly assign the first batch
        processed_df = head_chunk_result
    else:
        # Align the indices for addition with fill_value=0
        processed_df = processed_df.add(head_chunk_result, fill_value=0)

    buffer_df =  chunk[chunk["query"] != head_id]
    last_id = head_id

processed_df.reset_index().to_csv(args.output, sep="\t", header=is_first_q, index=False, mode='a')
