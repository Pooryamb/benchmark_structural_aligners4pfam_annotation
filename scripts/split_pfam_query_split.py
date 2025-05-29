import os
from fasta2dict import fasta2dict
from dict2fasta import dict2fasta

base_dir = "."
query_fasta = fasta2dict(f"{base_dir}/data/raw/dbs/pfam_split_query/pfam.fasta")

grouped_fasta = {}
for seed_id, seq in query_fasta.items():
    fam = seed_id.split("-")[-1]
    seq_dict = grouped_fasta.get(fam, {})
    seq_dict[seed_id] = seq
    grouped_fasta[fam] = seq_dict

out_dir = f"{base_dir}/data/raw/dbs/pfam_split_query/grp_by_family/"
os.makedirs(out_dir, exist_ok=True)

for fam, seq_dict in grouped_fasta.items():
    dict2fasta(seq_dict, f"{out_dir}/{fam}.fasta")
