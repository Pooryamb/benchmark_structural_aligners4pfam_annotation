import pandas as pd

input_file = "./data/pfam_clust_reps.tsv"
output_file= "./data/pfam_fl_clust_reps.tsv"


input_df = pd.read_csv(input_file, sep="\t", header=None)
input_df[0] = "AF-" + input_df[0].str.split("-", expand=True)[0] + "-F1-model_v4"
input_df = input_df.drop_duplicates()
input_df.to_csv(output_file, sep="\t", header=None, index=None)
