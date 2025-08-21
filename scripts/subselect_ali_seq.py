import os
import glob
import multiprocessing as mp
import polars as pl

script_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(script_dir, "..")


def subselect_ali_seq(pairs_path):
    pairs_filename = os.path.basename(pairs_path)
    pairing_tool = pairs_filename.split("_")[0]
    batch_info = pairs_filename.split("_")[2]
    tools = ["reseek", "fs"]

    for tool in tools:
        os.makedirs(f"{base_dir}/tmp/alis/split_pf_seq/balanced/{tool}", exist_ok=True)
    exh_tool = [x for x in tools if x != pairing_tool][0] # The name of the tool by which the comprehensive alignments were derived

    output_path_template = "{base_dir}/tmp/alis/split_pf_seq/balanced/{pairing_tool}/{exh_tool}_{batch_info}.tsv"

    if pairing_tool == "fs":
        src_paths = {"input_1": f"{base_dir}/tmp/alis/split_pf_seq/fs_pref_{batch_info}.tsv",
                     "input_2": f"{base_dir}/tmp/alis/split_pf_seq/reseek_exh_{batch_info}.tsv"}
    else:
        src_paths = {"input_1": f"{base_dir}/tmp/alis/split_pf_seq/reseek_sens_{batch_info}.tsv",
                     "input_2": f"{base_dir}/tmp/alis/split_pf_seq/fs_exh_{batch_info}.tsv"}

    src_paths["output_1"] = output_path_template.format(base_dir=base_dir, pairing_tool=pairing_tool, exh_tool=pairing_tool, batch_info=batch_info)
    src_paths["output_2"] = output_path_template.format(base_dir=base_dir, pairing_tool=pairing_tool, exh_tool=exh_tool, batch_info=batch_info)

    pairs_df = pl.read_csv(pairs_path, separator="\t", has_header=False)
    for index in [1, 2]:
        comp_df = pl.read_csv(src_paths[f"input_{index}"], separator="\t", has_header=False)
        sub_df = pairs_df.join(comp_df, on=["column_1", "column_2"], how="inner")
        sub_df.write_csv(src_paths[f"output_{index}"], separator="\t", include_header=False)


# Earilier, we made a list of pairs of structurally similar pairs with equal number of TPs and FPs
pair_paths = glob.glob(f"{base_dir}/tmp/structurally_similar_pairs/*.tsv")

#for pairs_path in pair_paths:
#    print(f"processing {pairs_path}")
#    subselect_ali_seq(pairs_path)

with mp.Pool(2) as pool:
    pool.map(subselect_ali_seq, pair_paths)
