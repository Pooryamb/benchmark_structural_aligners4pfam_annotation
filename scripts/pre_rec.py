import pandas as pd

def pre_rec(labeled_df):
    """calculates the precision and recall from a labeled dataframe.
    recall shows the fraction of total true positives captured.
    precision shows the fraction of predicted positives that are true positives."""
    tp = len(labeled_df[labeled_df['label'] == 1])
    fp = len(labeled_df[labeled_df['label'] == 0])

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / len(labeled_df) if len(labeled_df) > 0 else 0

    return precision, recall
