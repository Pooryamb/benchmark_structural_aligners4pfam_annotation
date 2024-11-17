#!/usr/bin/env python3

import glob
import os
import argparse

parser = argparse.ArgumentParser(description="Convert the tsv files of a database to a reseek cal database!")
# Add an optional argument for the file path
parser.add_argument(
    '--input_basename',
    type=str,
    required=True,  # Make this argument required
    help="""The path to the basename of the tsv files.
    It should not include the extension. Only the name should be here"""
    )
parser.add_argument(
        '--output_basename',
        type=str,
        help='The path to the basename of output'
    )

# Parse the arguments
args = parser.parse_args()

basename_path = args.input_basename

if args.output_basename is None:
    out_path = basename_path
else:
    out_path = args.output_basename

with open(f"{basename_path}.tsv") as fa_file, \
     open(f"{basename_path}_ca.tsv") as ca_file, \
     open(f"{basename_path}_h.tsv") as h_file, \
     open(f"{out_path}.cal", 'w') as ofile:
    for h_line in h_file:
        if not(h_line.strip()):
            continue
        ca_line = ca_file.readline()
        fa_line = fa_file.readline()
        ofile.write(f">{h_line.strip().split()[1]}\n")
        seq_list = list(fa_line.strip().split("\t")[1])
        seq_len = len(seq_list)
        coords = ca_line.strip().split("\t")[1].split(",")
        ca_x, ca_y, ca_z = coords[:seq_len], coords[seq_len:2*seq_len], coords[2*seq_len:3*seq_len]
        for i in range(seq_len):
            ofile.write(f"{seq_list[i]}\t{ca_x[i]}\t{ca_y[i]}\t{ca_z[i]}\n")
