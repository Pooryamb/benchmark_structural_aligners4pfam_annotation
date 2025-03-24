import pandas as pd
import numpy as np


def find_optimal_bins(df, stratifier_col, perf_func, initial_bins=20, final_bins=4):
    """
    Takes a dataframe as input, bins the rows based on the stratifier_col column. Binning is done based on
    quantiles. In the beginning, it uses the number of bins = initial_bins. It recursively decreases the
    number of bins by merging bins of similar performance until final_bins are left. After binning, the performance
    on each bin is calculated using the perf_func
    """
    df['range'] = pd.qcut(df[stratifier_col], q=initial_bins)  # Create 20 bins with ~equal frequencies

    # Step 2: Group by the quantile-based ranges and apply the `perf_func` perf_function
    group_stats = (
        df.groupby('range')
        .apply(perf_func)
        .reset_index()
        .rename(columns={0: 'performance'})
    )
    group_stats['range'] = group_stats['range'].astype(object)

    while len(group_stats) > final_bins:

        group_stats['diff'] = group_stats['performance'].diff().abs()
        min_diff_index = group_stats['diff'].idxmin()
        if min_diff_index > 0:
            left_range = group_stats.iloc[min_diff_index - 1]['range']
            right_range = group_stats.iloc[min_diff_index]['range']

            new_range = pd.Interval(left=left_range.left, right=right_range.right, closed='right')
            merged_data = df[df[stratifier_col].between(new_range.left, new_range.right)]
            new_perf = perf_func(merged_data)
            group_stats.loc[min_diff_index - 1, 'range'] = new_range
            group_stats.loc[min_diff_index - 1, 'performance'] = new_perf

            # Drop the current row (min_diff_index)
            group_stats = group_stats.drop(index=min_diff_index).reset_index(drop=True)

        group_stats['diff'] = group_stats['performance'].diff().abs()
    # Step 4: The final 4 ranges with their corresponding performance outputs
    final_ranges = group_stats.drop(columns=['diff'])
    return final_ranges
