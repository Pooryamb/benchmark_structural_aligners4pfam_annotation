def parse_stockholm_in_batches(filename):
    """
    Parses a Stockholm file and yields results for each entry (separated by "//").

    Args:
        filename (str): Path to the Stockholm file.

    Yields:
        dict: A dictionary with seed IDs as keys and aligned sequences as values for the current entry.
        str: The "AC" fields for the current entry.
    """
    current_seed_dict = {}  # Dictionary for seed IDs and aligned sequences
    current_ac_string = ""  # String for concatenated "AC" fields

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()  # Remove leading/trailing whitespace

            # Ignore empty lines or lines with comments that don't involve "#=GF AC"
            if line.startswith("#"):
                if line.startswith("#=GF AC"):  # Extract "AC" field
                    current_ac_string += line.split()[-1].strip()  # Extract and concatenate AC
                continue

            # Detect end of an entry
            if line == "//":
                # Yield the current entry's data
                yield current_seed_dict, current_ac_string
                # Reset for next entry
                current_seed_dict = {}
                current_ac_string = ""
                continue

            # Process a sequence line (seed ID and aligned sequence)
            parts = line.split()
            if len(parts) == 2:  # Ensure valid format (skip invalid lines)
                seed_id, aligned_seq = parts
                seed_id = seed_id.strip()
                aligned_seq = aligned_seq.strip().upper().replace(".", "-") #This is required for the correspondence mapping step
                # Add or append the aligned sequence to the dictionary
                current_seed_dict[seed_id] = aligned_seq
