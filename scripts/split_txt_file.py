# Generate a code snippet taking a tsv file as input and splitting it into files with less than N lines. N is provided as an argument.The default 
# value for N is 1 billion lines. The output will be stored in a directory provided as an argument. 
# Use linux split command to do the splitting.
import argparse
import os

parser = argparse.ArgumentParser(description="This script splits a large TSV file into smaller files with a specified maximum number of lines.")
parser.add_argument("--input", type=str, required=True, help="The path to the input TSV file")
parser.add_argument("--output_dir", type=str, help="The directory where the split files will be stored")
parser.add_argument("--max_lines", type=int, default=10_000_000, help="The maximum number of lines per split file (default: 10 million)")
args = parser.parse_args()

# If input ends with .tsv, make default output_dir based on input filename and store the split files in a directory named after the input file without extension
if not args.output_dir:
    base_name = os.path.basename(args.input)
    dir_name = os.path.splitext(base_name)[0]
    args.output_dir = os.path.join(os.path.dirname(args.input), dir_name)


os.makedirs(args.output_dir, exist_ok=True)

# Construct the split command
split_command = f"split -l {args.max_lines} --additional-suffix=.tsv {args.input} {os.path.join(args.output_dir, 'part_')}"
os.system(split_command)
