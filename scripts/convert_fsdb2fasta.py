# read the arguments from the command line using argparse
# I need to parse four arguments: db_path, subdb_type, output_path, foldseek_bin_path
import argparse
import os
import subprocess


parser = argparse.ArgumentParser(description="convert a foldseek databse to a fasta file")
parser.add_argument("--db_path", type=str, required=True, help="path to the foldseek database, the path should not include the subdb_type suffix")
parser.add_argument("--subdb_type", type=str, default="all", help="which subdatabase to convert. Options: 'all', '', '_ss', '_ca'. Default is 'all'")
parser.add_argument("--output_path", type=str, help="path to the basename of the output fasta, by default it will be named as <db_path>. The subdb_type and .fasta suffix will be added automatically")
parser.add_argument("--foldseek_bin_path", type=str, default="foldseek", help="path to the foldseek binary, by default it is assumed to be in the PATH")

args = parser.parse_args()
db_path = args.db_path
subdb_type = args.subdb_type
output_path = args.output_path
foldseek_bin_path = args.foldseek_bin_path

    # determine which subdb types to process
if subdb_type == "all":
    subdb_types = ["", "_ss", "_ca"]
else:
    subdb_types = [subdb_type]

for subdb in subdb_types:
    if subdb == "_ca":
        subdb_ext_preconv = "_ca_f64"
    else:
        subdb_ext_preconv = subdb

    db_file = f"{db_path}{subdb}"
    if not os.path.exists(db_file):
        print(f"Database file {db_file} does not exist, skipping...")
        continue

    if subdb != "":
        print(f"create header file for subdb type {subdb}")
        cmd = [foldseek_bin_path, "lndb", f"{db_path}_h", f"{db_path}{subdb_ext_preconv}_h"]
        subprocess.run(cmd, check=True)

    if subdb == "_ca":
        cmd = [foldseek_bin_path, "compressca", db_path, f"{db_path}{subdb_ext_preconv}", "--coord-store-mode", "3"]
        subprocess.run(cmd, check=True)

    if output_path:
        out_file = f"{output_path}{subdb}.fasta"
    else:
        out_file = f"{db_path}{subdb}.fasta"

    cmd = [foldseek_bin_path, "convert2fasta", f"{db_path}{subdb_ext_preconv}", out_file]
    subprocess.run(cmd, check=True)
    if subdb == "_ca":
        # remove the temporary databases created
        cmd = [foldseek_bin_path, "rmdb", f"{db_path}_ca_f64"]
        subprocess.run(cmd, check=True)
        cmd = [foldseek_bin_path, "rmdb", f"{db_path}_ca_f64_h"]
        subprocess.run(cmd, check=True)
    print(f"Converted {db_file} to {out_file}")
