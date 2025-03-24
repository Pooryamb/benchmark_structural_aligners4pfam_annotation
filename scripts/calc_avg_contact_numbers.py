from io import StringIO
import pickle
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import numpy as np

def parse_cal(path):
    content = open(path).read().lstrip(">").split("\n>")
    content_dict = {part.split("\n")[0].replace(".cif", ""): pd.read_csv(StringIO("\n".join(part.split("\n")[1:])), sep="\t", header=None)[[1, 2, 3]].rename(columns={1:0, 2:1, 3:2}) for part in content}
    return content_dict

def calculate_contact_numbers(df, cutoff=8.0):
    """
    Calculate the contact number for each residue in a protein.

    Parameters:
        df (pd.DataFrame): Input DataFrame with x, y, z columns representing coordinates.
        cutoff (float): Distance threshold for considering residues to be in contact.
    Returns:
        pd.DataFrame: Original DataFrame with an additional 'contact_number' column.
    """
    # Convert DataFrame to a NumPy array for fast computation
    coords = df.iloc[:, :3].to_numpy()  # Assumes x, y, z are in the first three columns
    num_residues = len(coords)  # Number of residues in the protein

    # Initialize an array to store contact numbers
    contact_numbers = np.zeros(num_residues)

    # Compute pairwise distances
    for i in range(num_residues):
        # Compute distances from residue `i` to all other residues
        distances = np.linalg.norm(coords[i] - coords, axis=1)

        # Count residues within the cutoff distance (exclude self-contact with -1)
        contact_numbers[i] = np.sum(distances < cutoff) - 1  # -1 to exclude residue i itself

    # Add contact numbers as a new column to the DataFrame
    return contact_numbers


input_path = "./data/raw/dbs/pfam_cif_cut_clust/pfam.cal"
in_dict_fmt = parse_cal(input_path)
with ProcessPoolExecutor() as executor:
    # Parallelize the dictionary comprehension
    in_dict_fmt = {key: result for key, result in zip(
        in_dict_fmt.keys(),
        executor.map(calculate_contact_numbers, in_dict_fmt.values())
    )}

in_dict_fmt = {x:y.mean() for x,y in in_dict_fmt.items()}

with open("data/processed/avg_contact_num.tsv", "w") as ofile:
    ofile.write("seed_id\tavg_contact_num\n")
    for key, value in in_dict_fmt.items():
        ofile.write(f"{key}\t{value}\n")
