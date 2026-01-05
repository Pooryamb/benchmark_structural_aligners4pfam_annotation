import pandas as pd
from overlap_ratio import overlap_ratio

def label_hits(ali_df, gs_df):
    """labels the hits in ali_df based on ground truth in gs_df.
    The ali_df should have columns:
    - unip_id
    - qstart
    - qend
    - pred_pfam
    The gs_df should have columns:
    - unip_id
    - start
    - end
    - db_id
    Returns a dataframe with an additional column 'label' where:
    - 1: true positive
    - 0: false positive
    - -1: new prediction
    """
    merged_df = ali_df.merge(gs_df, left_on=['unip_id'], right_on=['unip_id'], how='right', suffixes=('_ali', '_gs'))
    na_rows = merged_df["query"].isna()
    merged_df.loc[na_rows, "pred_fam_ali"] = "PF_____"
    merged_df.loc[na_rows, "qstart_ali"] = -2
    merged_df.loc[na_rows, "qend_ali"] = -1

    merged_df['overlap'] = merged_df.apply(lambda row: overlap_ratio([row['qstart_ali'], row['qend_ali']], [row['qstart_gs'], row['qend_gs']]), axis=1)
    merged_df["label"] = -1 # New predictions
    merged_df.loc[(merged_df['pred_fam_ali'] == merged_df['pred_fam_gs']) & (merged_df['overlap'] >= 0.25), 'label'] = 1 # True positive
    merged_df.loc[(merged_df['pred_fam_ali'] != merged_df['pred_fam_gs']) & (merged_df['overlap'] >= 0.25), 'label'] = 0 # False positive
    merged_df = merged_df.sort_values(by=['unip_id', 'qstart_gs', 'qend_gs', "label"], ascending=False)
    non_red_hits = merged_df.drop_duplicates(subset=['unip_id', 'qstart_gs', 'qend_gs'], keep='first')
    return non_red_hits.reset_index(drop=True)
