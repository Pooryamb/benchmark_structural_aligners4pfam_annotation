#!/usr/bin/env python3
import os
import argparse

parser = argparse.ArgumentParser(description="Convert a foldseek database to a reseek cal database")
# Add an optional argument for the file path
parser.add_argument(
    '--input_basename',
    type=str,
    required=True,  # Make this argument required
    help="""The path to the basename of the foldseek database
    It should not include the extension. Only the name should be here"""
    )

parser.add_argument(
        '--output_basename',
        type=str,
        help='The path to the basename of output'
    )

parser.add_argument(
        '--foldseek',
        type=str,
        help='The path to the executable foldseek program'
    )
    
# Parse the arguments
args = parser.parse_args()

inp_basename = args.input_basename

if args.output_basename is None:
    out_basename = inp_basename
else:
    out_basename = args.output_basename
    
if args.foldseek is None:
    foldseek = "foldseek"
else:
    foldseek = args.foldseek

os.system(f"{foldseek} lndb {inp_basename}_h {inp_basename}_ca_h")
os.system(f"{foldseek} convert2fasta {inp_basename} {out_basename}.fasta")
os.system(f"{foldseek} compressca {inp_basename} {inp_basename}_ca_f64 --coord-store-mode 3")
os.system(f"{foldseek} lndb {inp_basename}_h {inp_basename}_ca_f64_h")
os.system(f"{foldseek} convert2fasta {inp_basename}_ca_f64 {out_basename}_ca.fasta")
os.system(f"{foldseek} rmdb {inp_basename}_ca_f64")
os.system(f"{foldseek} rmdb {inp_basename}_ca_f64_h")


with open(f"{out_basename}_ca.fasta") as ca_file, \
     open(f"{out_basename}.fasta") as seq_file, \
     open(f"{out_basename}.cal", 'w') as ofile:
    
    while True:
    
        header_line = seq_file.readline()
        if not header_line.strip():
            break
        header = header_line.strip().strip(">")
        
        seq_list = list(seq_file.readline().strip())
        seq_len = len(seq_list)
        _ = ca_file.readline()
        ca_line = ca_file.readline()
        coords = ca_line.strip().split(",")
        
        ofile.write(f">{header}\n")    
        ca_x, ca_y, ca_z = coords[:seq_len], coords[seq_len:2*seq_len], coords[2*seq_len:3*seq_len]
        for i in range(seq_len):
            ofile.write(f"{seq_list[i]}\t{ca_x[i]}\t{ca_y[i]}\t{ca_z[i]}\n")
