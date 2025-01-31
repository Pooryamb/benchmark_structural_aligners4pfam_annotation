import pandas as pd

input_path  = "./data/processed/pfam_plddt.json"
output_path = "./data/processed/pfam_avg_plddt.tsv"

plddt_list = []

with open(input_path) as input_file:
    for line in input_file:
        if len(line) > 3:
            parts = line.split(" : ")
            seq_id = parts[0].lstrip(",").replace('"', '')
            exec("a = " + parts[1].strip())
            avg_plddt = round(sum(a)/len(a), 2)
            plddt_list.append([seq_id, avg_plddt])
plddt_df = pd.DataFrame(plddt_list)
plddt_df.to_csv(output_path, sep="\t", header=None, index=None)
