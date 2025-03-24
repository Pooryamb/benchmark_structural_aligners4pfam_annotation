import json
import glob
import os
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from tqdm import tqdm

from headers import fam_ali_headers
from find_pw_correspondence import find_pw_correspondence #Taking two sequences as input and the start of alignment on each sequnce(optional), it gives the mapping
from parse_stockholm import parse_stockholm

# load the pfam alignments as a dictionary
base_dir = "./"
pfam_dict, _ = parse_stockholm(open(f"{base_dir}/data/processed/aligned_pfam_clust.sto").read())
def unify_ali_dict(seq_dict):
    unified = {x:y.upper().replace(".", "-") for x,y in seq_dict.items()}
    return unified
pfam_dict = {x:unify_ali_dict(y) for x,y in pfam_dict.items()}

# load special residues of the seeds, such as pockets, conserved, etc.
special_residues_files = glob.glob(f"{base_dir}/tmp/residue_features/*.json") # It is assumed that the special residues whose alignment is of our interest is accessible in this directory
special_residues = {os.path.basename(x).replace(".json", ""): json.loads(open(x).read()) for x in special_residues_files}


tools = ["fs", "fs3di", "rs", "tm", "mm"]
header_line = "query\ttarget\tcorrect\tall\n"
def check_residue_alignment(family):
    """Takes the pfam family as input, and compares the residue level alignment for all members of the family with the pfam curated MSAs. It outputs multiple files, one for the overall
    residue alignment, others for residue alignment of special residues, such as conserved residues or residues located inside the pockets"""
    gs_msa = pfam_dict[family]
    for tool in tools:
        pwa = pd.read_csv(f"{base_dir}/tmp/alis/fam_alis/{tool}/{family}.tsv", sep="\t", header=None, names=fam_ali_headers[tool])
        pwa = pwa[pwa["query"]!=pwa["target"]]
        all_file = open(f"{base_dir}/tmp/intrafam_residue_alignment_counts/{tool}/all/{family}.tsv", 'w')
        special_files = {kind:open(f"{base_dir}/tmp/intrafam_residue_alignment_counts/{tool}/{kind}/{family}.tsv", 'w') for kind in special_residues.keys() }
        all_file.write(header_line)                                     # write header for the report of all residues
        _ = [f.write(header_line) for f in special_files.values()] # write header for the report of special residues
        for index, row in pwa.iterrows():
            pwa_mapping = find_pw_correspondence(row["qaln"], row["taln"], row["qstart"], row["tstart"])  # finds mapping between query and target for the selected tool
            gs_mapping = find_pw_correspondence(gs_msa[row["query"]], gs_msa[row["target"]])              # finds mapping between query and target based on Pfam MSAs
            correctly_mapped = set(pwa_mapping).intersection(set(gs_mapping))                             # finds the number of correctly mapped residues
            query, target = row["query"], row["target"]
            all_file.write(f"{query}\t{target}\t{len(correctly_mapped)}\t{len(gs_mapping)}\n")
            for kind, file in special_files.items():
                target_special_residues = special_residues[kind].get(target, [])
                special_coorectly_aligned = [x for x in correctly_mapped if x[1] in target_special_residues]
                file.write(f"{query}\t{target}\t{len(special_coorectly_aligned)}\t{len(target_special_residues)}\n")

    all_file.close()
    for file in special_files.values():
        file.close()

with ProcessPoolExecutor() as executor:
    # Parallelize the dictionary comprehension
    executor.map(check_residue_alignment, pfam_dict.keys())
