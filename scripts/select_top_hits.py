import os
from concurrent.futures import ProcessPoolExecutor
import glob
import pandas as pd
from headers import headers


def read_alns_in_chunks(path):
    """Takes the path to a file as input and yields chunks of data as output"""
    base_name = os.path.basename(path)
    tool = base_name.split("_")[0]
    if tool == "hmmscan":
        delim = r"\s+"
    else:
        delim = "\t"
    header = headers[tool]
    for subdf in pd.read_csv(path, sep=delim, header=None, names=header, comment="#", chunksize=10000):
        yield subdf



def select_top_hit(input_path, output_path = ""):
    """Selects the top hit for each query"""
    if not(output_path):
        output_path = input_path.replace("alis/split_pf", "first_hits")
    selected_rows = []
    for chunk in read_alns_in_chunks(input_path):
        selected_rows.append(chunk.drop_duplicates("query"))
    first_hits =  pd.concat(selected_rows).drop_duplicates("query")
    first_hits.to_csv(output_path, sep="\t", index=None)



if __name__ == "__main__":
    file_paths = [x for x in glob.glob("./tmp/alis/split_pf/*.tsv") if "_large_file" not in x]
    with ProcessPoolExecutor() as executor:
        executor.map(select_top_hit, file_paths)
