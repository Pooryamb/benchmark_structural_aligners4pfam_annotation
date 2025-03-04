def fasta2dict(path):
    with open(path) as ifile:
        content_parts = ifile.read().strip().lstrip(">").split("\n>")
        fasta_dict = {x.split("\n")[0]: "".join(x.split("\n")[1:]) for x in content_parts}
        return fasta_dict
