# Use argparse to parse 3 arguments: --fasta_path, --output_path, cut_coordinates
import argparse
argparser = argparse.ArgumentParser(description="takes the prefix path to fasta files corresponding to a foldseek database, cuts the sequences according to cut coordinates, and saves them to output path")
argparser.add_argument("--fasta_path", type=str, required=True, help="path to the fasta prefix files. It expects to find files with .fasta, _ss.fasta, and _ca.fasta suffixes")
argparser.add_argument("--output_path", type=str, required=True, help="path to save the cut sequences in tsv format. This one also should be a prefix path, as four files will be created with .tsv, _h.tsv, _ss.tsv, and _ca.tsv suffixes")
argparser.add_argument("--cut_recipe_path", type=str, help="path to a tsv file with four columns: sequence_id, family, start, end. The sequences will be cut according to these coordinates (1-based, inclusive)")
args = argparser.parse_args()

cut_recipe_path = args.cut_recipe_path
output_path = args.output_path
fasta_path = args.fasta_path


#!/usr/bin/env python3
# With special thanks to Dr. Milot Mirdita for providing the instructions on cutting a Foldseek database
def cut_normalseq(seq, domstart, domend):
    return seq[domstart-1:domend]

def cut_caseq(seq, domstart, domend):
    pdb_coords = seq.strip().split(",")
    res_count = len(pdb_coords)//3
    x_coords = pdb_coords[domstart - 1: domend]
    y_coords = pdb_coords[res_count + domstart-1: res_count + domend]
    z_coords = pdb_coords[2*res_count + domstart-1: 2*res_count + domend]
    all_coords = x_coords + y_coords + z_coords
    return ",".join(map(str, all_coords))

def extract_unipid_from_fasta(header):
    return header.split("-")[1]

cut_coords = {}
with open(cut_recipe_path) as seed_coords:
    for line in seed_coords:
        unipid, pfid, domstart, domend  = line.strip().split("\t")
        domstart, domend = int(domstart), int(domend)
        cut_positions = cut_coords.get(unipid, [])
        cut_positions.append([domstart, domend, pfid])
        cut_coords[unipid] = cut_positions

non_unipids= []
i=0
with open(f"{fasta_path}.fasta") as seq_file, \
     open(f"{fasta_path}_ss.fasta") as ss_seq_file, \
     open(f"{fasta_path}_ca.fasta") as ca_seq_file, \
     open(f"{output_path}.tsv", 'w') as cut_seq_file, \
     open(f"{output_path}_ss.tsv", 'w') as cut_ss_seq_file, \
     open(f"{output_path}_ca.tsv", 'w') as cut_ca_seq_file, \
     open(f"{output_path}_h.tsv", 'w') as cut_h_file:
    while True:
        header_line = seq_file.readline() # expects the header lines to be followed by sequence lines. Each sequence is in a single line
        if not header_line.startswith(">"):
            break
        else:
            unipid = extract_unipid_from_fasta(header_line)
        fl_seq = seq_file.readline().strip()
        _, __ = ss_seq_file.readline(),  ca_seq_file.readline()
        fl_ss_seq, fl_ca_seq = ss_seq_file.readline().strip(),  ca_seq_file.readline().strip()

        if unipid not in cut_coords:
            non_unipids.append(unipid)
            continue

        for cut_positions in cut_coords[unipid]:
            domstart, domend, pfid = cut_positions
            cut_seq_id = f"{unipid}_{domstart}_{domend}_{pfid}"
            cut_seq = cut_normalseq(fl_seq, domstart, domend)
            cut_ss_seq= cut_normalseq(fl_ss_seq, domstart, domend)
            cut_ca_seq=cut_caseq(fl_ca_seq, domstart, domend)

            cut_h_file.write(f"{i}\t{cut_seq_id}\n")
            cut_ca_seq_file.write(f"{i}\t{cut_ca_seq}\n")
            cut_ss_seq_file.write(f"{i}\t{cut_ss_seq}\n")
            cut_seq_file.write(f"{i}\t{cut_seq}\n")
            i+=1
