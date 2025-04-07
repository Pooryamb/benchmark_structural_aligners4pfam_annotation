import glob
import pandas as pd
import json


unip_data = pd.concat([pd.read_csv(x, sep="\t") for x in glob.glob("./data/raw/unip_sites_b*.tsv")])[["From", "Sequence"]]
unip_data = unip_data.rename(columns={"From": "unip_id", "Sequence": "seq"})
unip_data = unip_data.set_index("unip_id")

fl_active_sites = json.load(open("./data/processed/pfam_clust_fl_active_sites_unip_version.json"))
fl_binding_sites = json.load(open("./data/processed/pfam_clust_fl_binding_sites_unip_version.json"))

def find_all_occurrences(main_string, substring):
    """Finds all starting indices of a substring in a main string."""
    indices = []
    start_index = 0
    while True:
        index = main_string.find(substring, start_index)
        if index == -1:  # No more occurrences found
            break
        indices.append(index)
        start_index = index + 1 # Move start_index to search after the found occurrence
    return indices

def map_sites2seeds(fl_seq, fl_sites, dom_seq):
    """Takes the full-length sequence + its sites + domain sequence as input and returns the corresponding sites on the domain"""
    domain_offsets = find_all_occurrences(fl_seq, dom_seq)
    all_sites = set()
    for offset in domain_offsets:
        fl_sites_py_ind_offset = [x-offset for x in fl_sites]
        eligible_sites_first_occ = [x for x in fl_sites_py_ind_offset if x > 0 and x <= len(dom_seq)]
        all_sites = all_sites.union(set(eligible_sites_first_occ))
    return sorted(list(all_sites))

seeds_binding_sites = {}
seeds_active_sites = {}

with open("./data/raw/dbs/pfam_fs_cut_clust/pfam.tsv") as seq_file, \
     open("./data/raw/dbs/pfam_fs_cut_clust/pfam_h.tsv") as h_file:
    for seq_line in seq_file:
        if not (seq_line.strip()):
            break
        h_line = h_file.readline()
        dom_seq = seq_line.strip().split()[1]
        header = h_line.strip().split()[1]
        unip_id = header.split("-")[0]
        if unip_id in fl_binding_sites:
            fl_seq = unip_data.loc[unip_id]["seq"]
            mapped_sites = map_sites2seeds(fl_seq, fl_binding_sites[unip_id], dom_seq)
            if mapped_sites:
                seeds_binding_sites[header] = mapped_sites
        if unip_id in fl_active_sites:
            fl_seq = unip_data.loc[unip_id]["seq"]
            mapped_sites = map_sites2seeds(fl_seq, fl_active_sites[unip_id], dom_seq)
            if mapped_sites:
                seeds_active_sites[header] = mapped_sites


with open("./data/processed/binding_sites_on_seeds_just_unip_source.json", 'w') as ofile:
    json.dump(seeds_binding_sites, ofile, indent=4)

with open("./data/processed/active_sites_on_seeds_just_unip_source.json", 'w') as ofile:
    json.dump(seeds_active_sites, ofile, indent=4)
