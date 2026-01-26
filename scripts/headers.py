mm_suit_def_header = ['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits']
reseek_tbl = ['query', 'target', 'qstart', 'qend', 'qlen', 'tstart', 'tend', 'tlen', 'fident', 'evalue', 'aq']
tm_tbl = ['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'alntmscore', 'qtmscore', 'ttmscore']
hmmscan = ["target", "t_acc", "query", "q_acc", "evalue", "bits", "bias", "evalue_dom", "bits_dom", "bias_dom", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "desc"]
headers = {"mm":mm_suit_def_header, "fs_cut": mm_suit_def_header, "cif_cut": mm_suit_def_header, \
           "fs": mm_suit_def_header, "reseek": reseek_tbl, "tm": tm_tbl, "hmmscan": hmmscan, \
           "fs_cut_cif_cut": mm_suit_def_header, "cif_cut_fs_cut": mm_suit_def_header}

ali_col_additions = ["qaln", "taln"]
ali_lddt_col_additions = ["lddtfull", "qaln", "taln"]
fam_ali_headers = {"mm":mm_suit_def_header+ali_col_additions,
               "fs": mm_suit_def_header+ali_lddt_col_additions,
               "fs3di": mm_suit_def_header+ali_lddt_col_additions,
               "rs": reseek_tbl+ali_col_additions,
               "tm": tm_tbl+ali_lddt_col_additions}

