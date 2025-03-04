import glob
from concurrent.futures import ProcessPoolExecutor
import os
from dict2fasta import dict2fasta

def extract_ss_from_af_cif(path):
    with open(path) as file:
        read_line = False
        for line in file:
            if line.startswith("_ma_target_ref_db_details.seq_db_align_end"):
                prot_len = int(line.strip().split()[1])
                prot_ss_state = prot_len * ["C"]
            if line.startswith("_struct_conf.pdbx_end_PDB_ins_code"):
                read_line = True
                continue
            if read_line:
                if line.startswith("#"):
                    break
                start, end, ss_type_comp = int(line.split()[2]), int(line.split()[9]), line.split()[6]
                if ss_type_comp.startswith("HELX"):
                    ss_state = "H"
                elif ss_type_comp.startswith("STRN"):
                    ss_state = "E"
                else:
                    ss_state = "C"
                for res_num in range(start-1, end):
                    prot_ss_state[res_num] = ss_state
    return "".join(prot_ss_state)


if __name__ == "__main__":
    cif_files = glob.glob("./tmp/pfam_fl/*.cif")
    with ProcessPoolExecutor() as executor:
        ss_state_dict = {os.path.basename(path).split("-")[1] : executor.submit(extract_ss_from_af_cif, path).result() for path in cif_files}
    dict2fasta(ss_state_dict, "./data/processed/ss_pfam.fasta")

