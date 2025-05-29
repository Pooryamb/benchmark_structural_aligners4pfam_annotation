import json
import glob
import os
from concurrent.futures import ProcessPoolExecutor
import itertools
import pandas as pd

from headers import fam_ali_headers
from find_pw_correspondence import find_pw_correspondence
from parse_stockholm import parse_stockholm

# Load the pfam alignments as a dictionary
def load_pfam_data(base_dir):
    pfam_dict, _ = parse_stockholm(open(f"{base_dir}/data/processed/aligned_pfam_clust.sto").read())
    return {x: {k: v.upper().replace(".", "-") for k, v in y.items()} for x, y in pfam_dict.items()}

# Load special residues of the seeds, such as pockets, conserved, etc.
def load_special_residues(base_dir):
    special_residues_files = glob.glob(f"{base_dir}/tmp/residue_features/*.json")
    return {os.path.basename(x).replace(".json", ""): json.loads(open(x).read()) for x in special_residues_files}

# Main residue alignment checking function
def check_residue_alignment(args):
    family, tool, pfam_dict, special_residues = args

    """Compares the residue-level alignment for all members of the family with the pfam curated MSAs."""
    gs_msa = pfam_dict[family]

    pwa = pd.read_csv(f"{base_dir}/tmp/alis/fam_alis/{tool}/{family}.tsv", sep="\t", header=None, names=fam_ali_headers[tool])
    pwa = pwa[pwa["query"] != pwa["target"]]

    header_line = "query\ttarget\tcorrect\tall\n"
    all_file_path = f"{base_dir}/tmp/intrafam_residue_alignment_counts/{tool}/all/{family}.tsv"
    special_files_paths = {kind: f"{base_dir}/tmp/intrafam_residue_alignment_counts/{tool}/{kind}/{family}.tsv" for kind in special_residues.keys()}

    # Write to files
    with open(all_file_path, 'w') as all_file:
        all_file.write(header_line)

        special_file_handlers = {kind: open(file_path, 'w') for kind, file_path in special_files_paths.items()}
        for f in special_file_handlers.values():
            f.write(header_line)  # write header for special residues

        for index, row in pwa.iterrows():
            pwa_mapping = find_pw_correspondence(row["qaln"], row["taln"], row["qstart"], row["tstart"])
            gs_mapping = find_pw_correspondence(gs_msa[row["query"]], gs_msa[row["target"]])
            correctly_mapped = set(pwa_mapping).intersection(set(gs_mapping))
            query, target = row["query"], row["target"]
            all_file.write(f"{query}\t{target}\t{len(correctly_mapped)}\t{len(gs_mapping)}\n")
            for kind, file in special_file_handlers.items():
                seed_special_residues = special_residues[kind].get(query, [])
                special_correctly_aligned = [x for x in correctly_mapped if x[0] in seed_special_residues]
                file.write(f"{query}\t{target}\t{len(special_correctly_aligned)}\t{len(seed_special_residues)}\n")

        for f in special_file_handlers.values():
            f.close()

# Execution block
if __name__ == "__main__":
    base_dir = "./"

    pfam_dict = load_pfam_data(base_dir)
    special_residues = load_special_residues(base_dir)

    tools = ["fs", "fs3di", "rs", "tm", "mm"]

    all_fam_tool_combinations = sorted(list(itertools.product(pfam_dict.keys(), tools)))

    # Create arguments list to pass to the executor
    args = [(family, tool, pfam_dict, special_residues) for family, tool in all_fam_tool_combinations]

    with ProcessPoolExecutor() as executor:
        executor.map(check_residue_alignment, args)
