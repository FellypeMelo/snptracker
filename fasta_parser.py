"""
FASTA Parser Utility

Provides functions to read DNA sequences from standard FASTA files.
"""

def read_fasta(file_path: str) -> str:
    """
    Reads the first DNA sequence from a FASTA file.

    Args:
        file_path: Path to the FASTA file.

    Returns:
        str: The DNA sequence.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    sequence = []
    started = False
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if started:
                    # We've already read one sequence, return it
                    break
                started = True
            else:
                sequence.append(line)
                started = True # In case there's no header (technically invalid FASTA but robust)
                
    return "".join(sequence)
