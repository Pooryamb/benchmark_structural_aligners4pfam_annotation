import os
from concurrent.futures import ProcessPoolExecutor
import glob
import pandas as pd

def delete_file(file_path):
    """
    Deletes a single file. Returns the file path if successful, or an error message.
    """
    os.remove(file_path)


def delete_files_in_parallel(file_paths):
    """
    Deletes multiple files in parallel.
    """
    with ProcessPoolExecutor() as executor:
        executor.map(delete_file, file_paths)


clust_members = pd.read_csv("./data/raw/dbs/pfam_fs_cut_clust/pfam_h.tsv", sep="\t", header=None)
all_cifs = pd.DataFrame({1: glob.glob("./tmp/cifs/*")})
all_cifs[1] = all_cifs[1].str.replace("./tmp/cifs/", "").str.replace(".cif.gz", "")
non_clust_members = all_cifs.merge(clust_members, on=1, indicator=True, how="left")
non_clust_members = non_clust_members[non_clust_members["_merge"]=="left_only"][[1]]
non_clust_members[1] = "./tmp/cifs/" + non_clust_members[1] + ".cif.gz"
delete_files_in_parallel(list(non_clust_members[1]))
