import pandas as pd
import glob

cif_files = glob.glob("./tmp/cifs/*.cif")
with open("./tmp/extra/filepaths4prank.txt", 'w') as ofile:
    for path in cif_files:
        ofile.write(f"{path}\n")
