"""
SNPTracker - FASTA Parser Utility

Provides functions to read DNA sequences from standard FASTA files.

FASTA format expected:
    >header line (any text after '>')
    SEQUENCE (may span multiple lines; lines are concatenated)
    >next_header
    ...

Public API:
    read_fasta(file_path)         — returns the first sequence as a plain str
    read_all_sequences(file_path) — returns list[tuple[header, sequence]]

Both functions raise FileNotFoundError if the file does not exist.
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
    sequences = read_all_sequences(file_path)
    if not sequences:
        return ""
    return sequences[0][1]


def read_all_sequences(file_path: str) -> list[tuple[str, str]]:
    """
    Reads all DNA sequences from a FASTA file.

    Args:
        file_path: Path to the FASTA file.

    Returns:
        list[tuple[str, str]]: List of (header, sequence) tuples,
            one per entry in the file.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    results = []
    current_header = None
    current_sequence = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    results.append((current_header, "".join(current_sequence)))
                current_header = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line)

    if current_header is not None:
        results.append((current_header, "".join(current_sequence)))

    return results
