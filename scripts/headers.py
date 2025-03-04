mm_suit_def_header = ['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits']
reseek_tbl = ['query', 'target', 'qstart', 'qend', 'qlen', 'tstart', 'tend', 'tlen', 'fident', 'evalue', 'aq']
tm_tbl = ['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'alntmscore', 'qtmscore', 'ttmscore']
headers = {"mm":mm_suit_def_header, "fs_cut": mm_suit_def_header, "cif_cut": mm_suit_def_header, "reseek": reseek_tbl, "tm": tm_tbl}
