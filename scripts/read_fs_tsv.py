import pandas as pd

def read_fs_tsv(path):
    """Reads the foldseek tsv file with default columns"""
    fs_header = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits".split(",")
    fs_df = pd.read_csv(path, sep="\t", names=fs_header)
    fs_df["target"] = fs_df["target"].str.replace(".cif", "")
    fs_df["pred_fam"] = fs_df["target"].str.split("-", expand=True)[3]
    fs_df["unip_id"] = fs_df["query"].str.split("-", expand=True)[1]
    return fs_df
