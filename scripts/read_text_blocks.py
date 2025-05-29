def read_text_blocks(file_path, block_start="HMMER3/f", block_end="//"):
    """Generator that yields blocks from a large file. The words at the start and end of each block must be specified."""
    current_block = []
    in_block = False

    # Open the file for reading
    with open(file_path, 'r') as file:
        for line in file:
            line = line.rstrip('\n')  # Remove newline character

            # Check if the line indicates the start of a block
            if line.startswith(block_start):
                # If we were already in a block, yield the previous block
                if in_block:
                    yield current_block
                    current_block = []  # Reset for the next block

                # Start a new block
                in_block = True

            # If we are in a block, collect lines
            if in_block:
                current_block.append(line)

            # Check if the line indicates the end of a block
            if line.startswith(block_end):
                if in_block:
                    yield current_block  # Yield the current block
                    current_block = []  # Reset for the next block
                    in_block = False

    # Handle the case of the last block if the file does not end with //
    if in_block and current_block:
        yield current_block
