import re
import sys
import tarfile
import io
from io import StringIO
import gzip
import os
import argparse


def extract_plddt_cif(fhandle):
    """
    Convert a structure in mmCIF format to PDB format, and cuts the structure
    ------
    This function has been created by adjusting the code available in
    https://f1000research.com/articles/7-1961
    """

    in_section, read_atom = False, False
    label_pos = 0
    labels = {}
    all_plddt = []
    last_resnum = 0
    for line in fhandle:
        if line.startswith('loop_'):  # start of section
            in_section = True

        elif line.startswith('#'):  # end of section
            in_section = False
            read_atom = False

        elif in_section and line.startswith('_atom_site.'):  # ATOM/HETATM
            read_atom = True
            labels[line.strip()] = label_pos
            label_pos += 1

        elif read_atom and line.startswith(('ATOM', 'HETATM')):  # convert
            fields = re.findall(r'[^"\s]\S*|".+?"', line)  # find enclosed ''

            fid = labels.get('_atom_site.label_seq_id')
            resnum = int(fields[fid])

            bfac_col = labels.get('_atom_site.B_iso_or_equiv')
            res_plddt = float(fields[bfac_col])

            if resnum > last_resnum:
                all_plddt.append(res_plddt)
                last_resnum = resnum

    return all_plddt


def extract_plddt_from_tarredgz(input, output):
    with tarfile.open(input, 'r|') as tfile,\
        open(output, 'w') as outf:
        outf.write("{\n")
        for i,t in enumerate(tfile):
            gzFile = tfile.extractfile(t)
            f = gzip.open(gzFile, 'rt')
            seed_id = os.path.basename(t.get_info()['name'].replace(".cif.gz", ''))

            cif_lines = f.read().split("\n")
            lddt_vec = extract_plddt_cif(cif_lines)
            if i!=0:
                outf.write(",")
            outf.write("\"" + seed_id + "\" : " + str(lddt_vec) + "\n")
        outf.write("}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = """This script can extract 
    plddt from multiple tarred gzipped cif files. It stores the output in the json format""" )

    parser.add_argument("--input", required=True, type=str, help="Please specify the path to the input file")
    parser.add_argument("--output", required=True, type=str, help="Please specify the path to the json file to store the output file")
    args = parser.parse_args()
    extract_plddt_from_tarredgz(args.input, args.output)

