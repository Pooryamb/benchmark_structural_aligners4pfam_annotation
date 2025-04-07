import glob
import pandas as pd
import json


def extract_sites(site_info_from_unip, site_type="binding"):
    """Gets the information in the Binding/Active site columns of a row as input and
    returns the binding/active sites as a list"""
    if site_type == "binding":
        separator = "BINDING "
    elif site_type == "active":
        separator = "ACT_SITE "
    parts = site_info_from_unip.split(separator)[1:]
    parts = [x.split(";")[0] for x in parts]
    row_sites = []
    for part in parts:
        if ".." in part:
            range_start, range_end = int(part.split("..")[0]), int(part.split("..")[1])
            row_sites = row_sites + list(range(range_start, range_end + 1))
        else:
            row_sites = row_sites + [int(part)]
    return row_sites

sites_df = pd.concat([pd.read_csv(path, sep="\t") for path in glob.glob("./data/raw/unip_sites_b*.tsv")])

b_df = sites_df[~(sites_df["Binding site"].isna())]  # Rows with annotated binding sites
a_df = sites_df[~(sites_df["Active site"].isna())]   # Rows with annotated active sites

binding_sites_dict = {}
for index, row in b_df.iterrows():

    try:
        parsed_content = extract_sites(row["Binding site"], site_type="binding")
    except:
        print("The following row has an issue")
        print(row)
    binding_sites_dict[row["From"]] = parsed_content


active_sites_dict = {}
for index, row in a_df.iterrows():
    try:
        parsed_content = extract_sites(row["Active site"], site_type="active")
    except:
        print("The following row has an issue")
        print(row)
    active_sites_dict[row["From"]] = parsed_content



with open("./data/processed/pfam_clust_fl_binding_sites_unip_version.json", 'w') as ofile:
    json.dump(binding_sites_dict, ofile, indent=4)

with open("./data/processed/pfam_clust_fl_active_sites_unip_version.json", 'w') as ofile:
    json.dump(active_sites_dict, ofile, indent=4)
