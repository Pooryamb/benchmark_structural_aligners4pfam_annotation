from pathlib import Path
import tarfile
import gzip
import io
import os
import pandas as pd
from cut_cif import cut_cif

try:
    # script location -> parent of the script's directory (project root)
    base_dir = Path(__file__).resolve().parent.parent
except NameError:
    # running in a notebook or interactive shell where __file__ is not defined
    base_dir = Path.cwd().resolve().parent


pf_data = pd.read_csv(f"{base_dir}/data/tb_pfam.tsv", sep="\t", header=None)

dict_format = {}
for index, row in pf_data.iterrows():
    unip_id = row[0]
    unip_doms = dict_format.get(unip_id, [])
    unip_doms.append((row[1], row[2], row[3]))  # (domain_id, start, end)
    dict_format[unip_id] = unip_doms

tar_path = f"{base_dir}/data/raw/tb_pdb/UP000008524_185431_TRYB2_v4.tar"
tb_pfam_dir = f"{base_dir}/data/raw/tb_pfam"
os.makedirs(tb_pfam_dir, exist_ok=True)

with tarfile.open(tar_path, mode='r') as tar:    # 'r' since tar contains .cif.gz members
    for member in tar:
        if member.isreg() and member.name.endswith('.cif.gz'):
            with tar.extractfile(member) as raw:    # raw is a file-like for the compressed bytes
                if raw is None:
                    continue
                # Wrap the gzip stream in a text wrapper for line-by-line reading
                with gzip.GzipFile(fileobj=raw) as gz:
                    with io.TextIOWrapper(gz, encoding='utf-8') as fh:
                        cif_lines = fh.readlines()
                        unip_id = member.name.split('-')[1]
                        if unip_id in dict_format:
                            for domain_info in dict_format[unip_id]:
                                domain_id, start, end = domain_info
                                cut_cif_content = cut_cif(cif_lines, start, end)
                                with open(f"{tb_pfam_dir}/{unip_id}-{start}-{end}-{domain_id}.cif", "w") as out_fh:
                                    out_fh.write(cut_cif_content)
