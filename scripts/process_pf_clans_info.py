import pandas as pd


df = pd.read_csv("./tmp/Pfam-A.clans.tsv", sep="\t", header=None)
df = df[[0,1]]
df.loc[df[1].isna(),1] = df.loc[df[1].isna(),0]
df = df.rename(columns={0:"pfam", 1: "clan"})
df.to_csv("./data/pfam_clan_info.tsv", sep="\t", index=None)

