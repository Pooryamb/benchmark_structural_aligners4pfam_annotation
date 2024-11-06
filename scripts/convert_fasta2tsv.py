import argparse
import os



parser = argparse.ArgumentParser(description="""This script takes the basename to a fasta file made from  a foldseek database as input and makes tsv
                                 files of the same database. .fasta, _ca.fasta, and _ss.fasta files should be available for the provided basename""")

parser.add_argument("--input_basename", type=str, required=True, help="The path to the basename of the fasta files.")
args = parser.parse_args()
basename = args.input_basename

with open(f"{basename}.fasta") as seq_file, \
    open(f"{basename}_ss.fasta") as ss_file, \
    open(f"{basename}_ca.fasta") as ca_file, \
    open(f"{basename}.tsv", 'w') as seq_tsv, \
    open(f"{basename}_ss.tsv", 'w') as ss_tsv, \
    open(f"{basename}_ca.tsv", 'w') as ca_tsv, \
    open(f"{basename}_h.tsv", 'w') as h_tsv :
    i = 0
    while True:

        h_line = seq_file.readline().strip(">")
        if not h_line:
            break
        seq_line = seq_file.readline()
        _, ss_line = ss_file.readline(), ss_file.readline()
        _, ca_line = ca_file.readline(), ca_file.readline()
        seq_tsv.write(f"{i}\t{seq_line}")
        h_tsv.write(f"{i}\t{h_line}")
        ca_tsv.write(f"{i}\t{ca_line}")
        ss_tsv.write(f"{i}\t{ss_line}")
