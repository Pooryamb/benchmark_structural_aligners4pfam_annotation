def dict2fasta(seq_dict, out_path):
    with open(out_path, 'w') as ofile:
        for seq_id, ss_seq in seq_dict.items():
            ofile.write(f">{seq_id}\n{ss_seq}\n")

