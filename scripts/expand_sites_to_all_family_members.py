import glob
import json
import pandas as pd

from parse_stockholm import parse_stockholm
from map_residues import map_poses_from_aligned_seqs2entries


parsed_msa, _ = parse_stockholm(open("./data/processed/aligned_pfam_clust.sto").read())

def convert_seed_data2fam_data(input_dict):
    """If a dictionary has information for seeds, it will convert it to family information, so that
    the keys would be family names, and each value will be a dictionary with keys of seed ids and original values"""
    fam_data_dict = {}
    for seed_id, seed_values in input_dict.items():
        fam_id = seed_id.split("-")[3]
        fam_specific_data = fam_data_dict.get(fam_id, {})
        fam_specific_data[seed_id] = seed_values
        fam_data_dict[fam_id] = fam_specific_data
    return fam_data_dict

gap_char = "."
def map_seed_poses2msa(msa_dict, all_sites_dict):
    """The function takes two the MSAs parsed by parse_sto and also the dictionary of all
    conserved columns as input and returns the the conserved columns in each family as output"""
    all_sites = {}
    # all_sites_data = {} # This was to check how often important sites align with each other
    for pf, fam_sites_dict in all_sites_dict.items():
        # fam_annotated_sites_num = [len(y) for y in fam_sites_dict.values()], # This was to check how often important sites align with each other
        # max_ann_col_count = sum(fam_annotated_sites_num)   # It happens when there is no overlap at all among the sites, # This was to check how often important sites align with each other
        # min_ann_col_count = max(fam_annotated_sites_num)   # This is the minimum number of possible important sites in msa, # This was to check how often important sites align with each other
        fam_sites = set()
        fam_seq_dict = msa_dict[pf]

        for seed_id, seed_sites in fam_sites_dict.items():
            aligned_seed_seq = fam_seq_dict[seed_id]
            seed_sites = fam_sites_dict[seed_id]
            seed_pointer = -1
            for msa_pointer in range(len(aligned_seed_seq)):
                if aligned_seed_seq[msa_pointer] != gap_char:
                    seed_pointer += 1
                    ind_1_seed_pointer = seed_pointer + 1
                    if ind_1_seed_pointer in seed_sites:
                        fam_sites.add(msa_pointer + 1) # The stored indices start from 1
        all_sites[pf] = sorted(list(fam_sites))
        # all_sites_data[pf] = [min_ann_col_count, max_ann_col_count, len(fam_sites)] # This was to check how often important sites align with each other
    return all_sites #, all_sites_data # This was to check how often important sites align with each other

def map_msa_cols2seed_cols(msa_dict, msa_cols, msa_cols_mapped2seeds):
    """The function takes two inputs, the MSA dictionary and the conserved columns of each family as input and
    returns a dictionary containing the important sites of each seed. The MSA dictionary should have been parsed
    using parse_sto function"""
    all_seeds_sites_dict = {}
    for pf, pf_sites in msa_cols.items():
        seq_dict = parsed_msa[pf]
        for seed_id in seq_dict:

            seed_sites = []
            for site in pf_sites:
                mapped_site = msa_cols_mapped2seeds[pf][seed_id].get(site, -1)
                #if seed_id == "P10539-137-328-PF02774":
                if mapped_site != -1:
                    seed_sites.append(mapped_site)

            all_seeds_sites_dict[seed_id] = seed_sites
    return all_seeds_sites_dict

active_sites_from_unip = json.loads(open("./data/processed/active_sites_on_seeds_just_unip_source.json").read())
fam_seed_active_sites = convert_seed_data2fam_data(active_sites_from_unip)
fam_active_sites = map_seed_poses2msa(parsed_msa, fam_seed_active_sites)

binding_sites_from_unip = json.loads(open("./data/processed/binding_sites_on_seeds_just_unip_source.json").read())
fam_seed_binding_sites = convert_seed_data2fam_data(binding_sites_from_unip)
fam_binding_sites = map_seed_poses2msa(parsed_msa, fam_seed_binding_sites)

msa_cols_mapped2seeds = {x : map_poses_from_aligned_seqs2entries(y) for x,y in parsed_msa.items()}


active_sites_all_seeds = map_msa_cols2seed_cols(parsed_msa, fam_active_sites, msa_cols_mapped2seeds)
with open("./tmp/residue_features/active_site.json", 'w') as ofile:
    json.dump(active_sites_all_seeds, ofile, indent=4)

binding_sites_all_seeds = map_msa_cols2seed_cols(parsed_msa, fam_binding_sites, msa_cols_mapped2seeds)
with open("./tmp/residue_features/binding_site.json", 'w') as ofile:
    json.dump(binding_sites_all_seeds, ofile, indent=4)
