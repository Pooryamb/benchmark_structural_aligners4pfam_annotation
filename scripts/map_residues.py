def map_poses_from_aligned_seqs2entries(seq_dict, gap_char="."):
    """Takes a dictionary containing aligned sequences as input and outputs the mapping between stockholm columns
    and residues of each member of dictionary"""
    mapping = {}
    for seed_id, seed_seq in seq_dict.items():
        mapping[seed_id] = {}
        seed_pointer = 0
        for sto_pointer in range(len(seed_seq)):
            if seed_seq[sto_pointer].upper() != gap_char:
                seed_pointer += 1
                mapping[seed_id][sto_pointer+1] = seed_pointer  # indices start from 1
    return mapping
