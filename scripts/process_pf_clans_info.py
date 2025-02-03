import pandas as pd


df = pd.read_csv("./data/raw/Pfam-A.clans.tsv", sep="\t", header=None)
df = df[[0,1]]
df.loc[df[1].isna(),1] = df.loc[df[1].isna(),0]
df = df.rename(columns={0:"pfam", 1: "clan"})
df.to_csv("./data/raw/pfam_clan_info.tsv", sep="\t", index=None)

