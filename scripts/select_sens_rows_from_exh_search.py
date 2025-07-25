import glob
import polars as pl



def select_exh_subset(exh_path):    
    sens_path = exh_path.replace("split_pf_seq", "split_pf").replace("exh", "sens")
    out_path =  exh_path.replace("exh", "sens")
    sens_df = pl.read_csv(sens_path, separator="\t",has_header=False, columns=[0, 1])
    
    exh_df = pl.read_csv(exh_path, separator="\t",has_header=False)
    sel_df = sens_df.join(exh_df, on=["column_1", "column_2"], how="inner")
    sel_df.write_csv(out_path, separator="\t", include_header=False)


exh_paths = glob.glob("./tmp/alis/split_pf_seq/reseek_exh_*.tsv")
for exh_path in exh_paths:
    select_exh_subset(exh_path)
