def overlap_ratio(interval1, interval2):
    """
    Calculate the overlap length divided by the length of the second interval.

    Args:
        interval1: tuple or list [start1, end1]
        interval2: tuple or list [start2, end2]

    Returns:
        float: overlap_length / length_of_interval2
    """
    start1, end1 = interval1
    start2, end2 = interval2

    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap_length = max(0, overlap_end - overlap_start)

    length2 = end2 - start2
    if length2 == 0:
        return 0.0  # Avoid division by zero

    return overlap_length / length2
