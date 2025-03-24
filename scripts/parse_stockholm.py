def parse_stockholm(sto_content):
    """takes a string that is in stockholm format and returns two dictionaries. The keys of the first dictionary are
    family accessions of Pfam members and its values are dictionaries of sequences. The keys of each sequence dictionary are
    the headers of the sequences and their values are the sequences. In this case, the aligned sequences are stored."""
    pfam_fam_seq_dict = {}
    acc2id = {}
    lines = sto_content.split("\n")
    for line in lines:
        if not(line.strip()):
            break
        elif line.startswith("# STOCKHOLM 1.0"):
            fam_seqs = {}
        elif line.startswith("#=GF ID   "):
            family_id = line.split()[-1]
        elif line.startswith("#=GF AC   "):
            family_acc = line.strip().split()[-1].split(".")[0]
            acc2id[family_acc] = family_id
        elif line[:2].isalnum():
            header = line.strip().split()[0]
            seq = line.strip().split()[1]
            fam_seqs[header] = seq
        elif line.startswith("//"):
            pfam_fam_seq_dict[family_acc] = fam_seqs
    return pfam_fam_seq_dict, acc2id


def write_parsed_sto(parsed_sto, acc2id, output_path):
    """Takes a parsed stockholm format, accession to id mapping, and an output path as input
    and writes all into a file whose path is specified by output_path"""
    with open(output_path, 'w') as ofile:
        for fam, fam_seq_dict in parsed_sto.items():
            ofile.write("# STOCKHOLM 1.0\n")
            ofile.write(f"#=GF ID   {acc2id[fam]}\n")
            ofile.write(f"#=GF AC   {fam}\n")
            for header, seq in fam_seq_dict.items():
                ofile.write(f"{header:<40}{seq}\n")
            ofile.write("//\n")


