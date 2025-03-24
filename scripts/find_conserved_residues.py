from collections import Counter
import math
import json
import pandas as pd
from parse_stockholm import parse_stockholm
from map_residues import map_poses_from_aligned_seqs2entries

base_dir = "."
parsed_sto, _ = parse_stockholm(open(f"{base_dir}/data/processed/aligned_pfam_clust.sto").read())
parsed_sto_populated = {x:y for x,y in parsed_sto.items() if len(y)>=3} # We only select the Pfam families with more than 3 members

def calculate_shannon_entropy(amino_acids, gap_char='.'):
    """
    Calculate the Shannon entropy for an amino acid sequence, considering or ignoring gaps.

    Parameters:
        amino_acids (list or str): A list or string of amino acids (residues) that may contain gaps.
        gap_char (str): The character representing gaps (default: "-").

    Returns:
        float: The Shannon entropy of the residues.
    """
    # If the majority of a column is made of gaps, the column will be ignored.
    non_gap_amino_acids = [residue.upper() for residue in amino_acids if residue != gap_char]
    if len(non_gap_amino_acids) < 3 or len(non_gap_amino_acids)/len(amino_acids) < 0.8:
        return None
    # Count the frequency of each amino acid (or gap if included)
    counts = Counter(non_gap_amino_acids)
    total = sum(counts.values())
    # Handle the case where the column is all gaps (e.g., after filtering)

    # Calculate relative probabilities and entropy
    entropy = 0
    for residue, count in counts.items():
        probability = count / total
        entropy -= probability * math.log2(probability)
    return entropy


# We look for the positions in aligned sequences where the Shannon entropy is below 0.5
ent_threshold = 0.5
coserved_positions = {}
for fam_acc, fam_dict in parsed_sto_populated.items():
    fam_ali_len = len(list(fam_dict.values())[0])

    map2seeds = map_poses_from_aligned_seqs2entries(fam_dict) # The mapping between aligned sequences and seeds are found
    for i in range(fam_ali_len):
        residues_i = [x[i] for x in fam_dict.values()]
        entropy_i = calculate_shannon_entropy(residues_i) # Entropy of the ith column is calculated if it is below the threshold, the corresponding location on seeds will be stores
        if not(entropy_i is None) and entropy_i < ent_threshold:
            for seed_id in fam_dict.keys():
                seed_conserved_positions = coserved_positions.get(seed_id, [])
                if not(map2seeds[seed_id].get(i+1, None) is None):
                    seed_conserved_positions.append(map2seeds[seed_id][i+1])
                coserved_positions[seed_id] = seed_conserved_positions
with open(f"{base_dir}/tmp/residue_features/conserved.json", 'w') as file:
    json.dump(coserved_positions, file, indent=4)  #The location of conserved positions is written
