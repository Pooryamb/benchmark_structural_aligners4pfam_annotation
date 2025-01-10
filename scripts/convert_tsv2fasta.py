import os
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("--input_basename", type=str, required=True, help="""The path to the basename of the tsv file. tsv files for header, 
amino acid sequence, ss sequence, and carbon alpha should be present""")
parser.add_argument("--output_basename", type=str, help="The path to the basename of the output. If not provided, it will create the fasta file in the same path as the input")

args = parser.parse_args()

basename_path = args.input_basename

if args.output_basename is None:
    out_path = basename_path
else:
    out_path = args.output_basename


with open(f"{out_path}.fasta", 'w') as ofile, \
     open(f"{basename_path}.tsv") as seq_file, \
     open(f"{basename_path}_h.tsv") as h_file:
    while True:
        seq_line = seq_file.readline()
        h_line = h_file.readline()
        if seq_line.strip() == "":
            break
        ofile.write(f">{h_line.strip().split()[1]}\n{seq_line.strip().split()[1]}\n")
