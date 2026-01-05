import argparse
argparser = argparse.ArgumentParser(description="takes the prefix path to tsv files corresponding to a foldseek database, and makes a foldseek database out of them")
argparser.add_argument("--tsv_path", type=str, required=True, help="path to the tsv prefix files. It expects to find files with .tsv, _ss.tsv, _h.tsv, and _ca.tsv suffixes")
argparser.add_argument("--output_path", type=str, help="path to save the foldseek database. If it is not provided, it will be saved in the same directory as the input tsv files with the same name (without suffixes)")
argparser.add_argument("--foldseek", default="foldseek", type=str, help="path to foldseek executable, default assumes it is in PATH")
args = argparser.parse_args()

tsv_path = args.tsv_path
output_path = args.output_path if args.output_path else tsv_path
foldseek = args.foldseek

import os

os.system(f"{foldseek} tsv2db {tsv_path}.tsv {output_path} --output-dbtype 0")
os.system(f"{foldseek} tsv2db {tsv_path}_ss.tsv {output_path}_ss --output-dbtype 0")
os.system(f"{foldseek} tsv2db {tsv_path}_h.tsv {output_path}_h --output-dbtype 12")
os.system(f"{foldseek} tsv2db {tsv_path}_ca.tsv {output_path}_ca --output-dbtype 12")
os.system(f"{foldseek} compressca {output_path} {output_path}_ca2 --coord-store-mode 2")
os.system(f"{foldseek} rmdb {output_path}_ca")
os.system(f"{foldseek} mvdb {output_path}_ca2 {output_path}_ca")
