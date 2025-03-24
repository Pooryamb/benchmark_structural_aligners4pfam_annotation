def find_pw_correspondence(seq1, seq2, seq1_start=1, seq2_start=1):
    """
    Find the pairwise correspondence between two sequences (with gaps). It assumes that read sequences are in upper case
    format and gaps are shown with '-' signs.

    Args:
        seq1: The first sequence (string) with residues and gaps ('-' or '.').
        seq2: The second sequence (string) with residues and gaps ('-' or '.').

    Returns:
        A list of tuples. Each tuple contains:
            (position in seq1, position in seq2)
    """

    correspondence = []  # List to store the correspondences
    pos1, pos2 = seq1_start, seq2_start    # Positions in seq1 and seq2 (respective indices WITHOUT gaps)

    # Traverse both sequences together
    for char1, char2 in zip(seq1, seq2):
        # If both characters are residues, record their correspondence
        if char1 not in '-.' and char2 not in '-.':
            correspondence.append((pos1, pos2))
        if char1 != '-' and char1 != '.':  # If seq1 has a residue (not a gap)
            pos1 += 1  # Increment position in seq1 (1-based indexing)
        if char2 != '-' and char2 != '.':  # If seq2 has a residue (not a gap)
            pos2 += 1  # Increment position in seq2 (1-based indexing)
    return correspondence
