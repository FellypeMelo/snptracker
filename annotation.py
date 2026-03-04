"""
SNPTracker - SNP Annotation Module

Provides codon translation and functional annotation of SNPs
using the standard genetic code.

Assumes reading frame starts at position 1 (frame +1).
"""

CODON_TABLE: dict[str, str] = {
    # Phe
    "TTT": "Phe", "TTC": "Phe",
    # Leu
    "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    # Ile
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile",
    # Met (start codon)
    "ATG": "Met",
    # Val
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    # Ser
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "AGT": "Ser", "AGC": "Ser",
    # Pro
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    # Thr
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    # Ala
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    # Tyr
    "TAT": "Tyr", "TAC": "Tyr",
    # Stop codons
    "TAA": "STOP", "TAG": "STOP", "TGA": "STOP",
    # His
    "CAT": "His", "CAC": "His",
    # Gln
    "CAA": "Gln", "CAG": "Gln",
    # Asn
    "AAT": "Asn", "AAC": "Asn",
    # Lys
    "AAA": "Lys", "AAG": "Lys",
    # Asp
    "GAT": "Asp", "GAC": "Asp",
    # Glu
    "GAA": "Glu", "GAG": "Glu",
    # Cys
    "TGT": "Cys", "TGC": "Cys",
    # Trp
    "TGG": "Trp",
    # Arg
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGA": "Arg", "AGG": "Arg",
    # Gly
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

_VALID_BASES = frozenset("ACGT")


def get_codon(sequence: str, position: int) -> str:
    """Returns the codon (triplet) that contains the given 1-indexed position.

    Assumes reading frame starts at position 1. Returns an empty string if
    the position falls within an incomplete (trailing) codon.

    Args:
        sequence: Full DNA sequence string (uppercase).
        position: 1-indexed SNP position within the sequence.

    Returns:
        str: 3-character codon, or "" if the codon is incomplete.

    Raises:
        ValueError: If position is <= 0 or exceeds the sequence length.
    """
    if position <= 0 or position > len(sequence):
        raise ValueError(
            f"Position {position} is out of range for sequence of "
            f"length {len(sequence)}."
        )

    if len(sequence) < 3:
        return ""

    codon_index = (position - 1) // 3
    start = codon_index * 3
    codon = sequence[start:start + 3]

    if len(codon) < 3:
        return ""

    return codon


def translate_codon(codon: str) -> str:
    """Translates a codon to its amino acid using the standard genetic code.

    Args:
        codon: 3-character DNA codon string (uppercase, bases A/C/G/T only).

    Returns:
        str: Amino acid 3-letter abbreviation (e.g., "Ala") or "STOP".

    Raises:
        ValueError: If codon is not a valid 3-character uppercase ACGT string.
    """
    if len(codon) != 3 or not all(b in _VALID_BASES for b in codon):
        raise ValueError(
            f"Invalid codon '{codon}'. Expected 3 uppercase bases (A/C/G/T)."
        )
    return CODON_TABLE[codon]


def annotate_snp(position: int, ref_sequence: str, alt_sequence: str) -> str:
    """Annotates the functional effect of a SNP on the encoded amino acid.

    Identifies the codon containing the SNP in both reference and alternate
    sequences, translates both codons, and classifies the change.

    Assumes reading frame starts at position 1. Returns "NON_CODING" when
    the position falls in an incomplete codon or sequences are too short.

    Args:
        position: 1-indexed position of the SNP.
        ref_sequence: Full reference DNA sequence (uppercase).
        alt_sequence: Full alternate (sample) DNA sequence (uppercase).

    Returns:
        str: One of "SYNONYMOUS", "NON_SYNONYMOUS", "NONSENSE", "NON_CODING".
    """
    ref_codon = get_codon(ref_sequence, position)
    alt_codon = get_codon(alt_sequence, position)

    if not ref_codon or not alt_codon:
        return "NON_CODING"

    try:
        ref_aa = translate_codon(ref_codon)
        alt_aa = translate_codon(alt_codon)
    except (ValueError, KeyError):
        return "NON_CODING"

    if alt_aa == "STOP":
        return "NONSENSE"
    if ref_aa == alt_aa:
        return "SYNONYMOUS"
    return "NON_SYNONYMOUS"
