import os
from parse_stockholm import parse_stockholm, write_parsed_sto

base_dir = "."
target_sto, acc2id = parse_stockholm(open(f"{base_dir}/data/raw/dbs/pfam_split_target/pfam.sto").read())

out_dir = f"{base_dir}/data/raw/dbs/pfam_split_target/grp_by_family/"
os.makedirs(out_dir, exist_ok=True)
for pf, seq_dict in target_sto.items():
    fam_sto = {pf:seq_dict}
    write_parsed_sto(fam_sto, acc2id, f"{out_dir}/{pf}.sto")
