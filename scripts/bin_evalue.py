import os
import argparse
import pandas as pd
import numpy as np



# Create the parser
parser = argparse.ArgumentParser(description="Finds the number of rows with each e-value")  

# Add arguments  
parser.add_argument('--input_path', type=str, help='The path to the input tsv file', required=True)
parser.add_argument('--output_path', type=str, help='The path to the processed bin files')

# Parse the arguments
args = parser.parse_args()
if args.output_path is None:
    file_name_base = os.path.basename(args.input_path)
    args.output_path = f"./data/processed/evalue_bins/{file_name_base}"

# Initialize an empty Series to store bin counts
global_bin_counts = pd.Series(dtype=int)

# Define a function to process each chunk
def process_chunk(chunk):
    global global_bin_counts

    # Select the second-to-last column
    col_name = chunk.columns[-2]

    # Apply the transformation round(-log(value))
    transformed = chunk[col_name].apply(lambda x: -round(np.log(x)) if x > 0 else 1000)

    # Count occurrences in each bin for the current chunk
    counts = transformed.value_counts()

    # Aggregate the counts into the global bin counts
    global_bin_counts = global_bin_counts.add(counts, fill_value=0)

# Read the TSV file in chunks
chunk_size = 100_000

if __name__ == "__main__":

    for chunk in pd.read_csv(args.input_path, sep='\t', chunksize=chunk_size):
        process_chunk(chunk)
# Convert the global_bin_counts to a sorted DataFrame for better readability
    result = global_bin_counts.sort_index().reset_index()
    result.columns = ['-round(log(evalue))', 'Count']
    result.to_csv(args.output_path, sep="\t", index=None)
