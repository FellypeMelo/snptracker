"""
SNPTracker - SNP Annotation Module

Provides codon translation and functional annotation of SNPs
using the standard genetic code (all 64 codons).

Assumptions:
    - Sequences are uppercase strings containing only A, C, G, T.
    - Position is always 1-indexed (forward strand).

Supported reading frames:
    Forward strand: +1 (default), +2, +3
    Reverse complement strand: -1, -2, -3

Annotation categories returned by annotate_snp():
    SYNONYMOUS     — substitution that does not change the amino acid
    NON_SYNONYMOUS — substitution that changes the amino acid
    NONSENSE       — substitution that introduces a stop codon (TAA/TAG/TGA)
    NON_CODING     — SNP falls in an incomplete codon, before frame start,
                     or sequence too short

Public API:
    CODON_TABLE              — dict mapping 3-mer → amino acid or "STOP"
    reverse_complement(seq)  — returns the reverse complement of a sequence
    get_codon(seq, pos, frame=1)
        — returns the triplet containing 1-indexed pos under the given frame
    translate_codon(codon)   — returns amino acid abbreviation or "STOP"
    annotate_snp(pos, ref, alt, frame=1)
        — annotates a SNP under the given reading frame
    is_in_cds(pos, cds_regions)
        — returns True if position is within any (start, end) CDS region
    annotate_snp_with_regions(pos, ref, alt, cds_regions=None, frame=1)
        — like annotate_snp but returns NON_CODING for positions outside CDS;
          cds_regions=None falls back to annotate_snp (backwards-compatible)
"""

_VALID_BASES = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    """Returns the reverse complement of a DNA sequence.

    Args:
        sequence: Uppercase DNA string (A/C/G/T only).

    Returns:
        str: Reverse complement sequence.

    Raises:
        ValueError: If sequence contains non-ACGT characters.
    """
    if any(b not in _VALID_BASES for b in sequence):
        raise ValueError(
            f"Sequence contains non-ACGT characters: '{sequence}'."
        )
    return sequence.translate(_COMPLEMENT)[::-1]


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


def get_codon(sequence: str, position: int, frame: int = 1) -> str:
    """Returns the codon (triplet) that contains the given 1-indexed position.

    For forward frames (+1, +2, +3), the reading frame starts at the given
    offset within the forward strand. For reverse frames (-1, -2, -3), the
    codon is looked up on the reverse complement strand.

    Args:
        sequence: Full DNA sequence string (uppercase).
        position: 1-indexed SNP position within the sequence.
        frame: Reading frame. One of {1, 2, 3, -1, -2, -3}. Default 1.

    Returns:
        str: 3-character codon, or "" if the position is NON_CODING in this
             frame (incomplete codon or before frame start).

    Raises:
        ValueError: If position is out of range or frame is invalid.
    """
    _VALID_FRAMES = {1, 2, 3, -1, -2, -3}
    if frame not in _VALID_FRAMES:
        raise ValueError(
            f"Invalid frame '{frame}'. Must be one of {_VALID_FRAMES}."
        )

    if position <= 0 or position > len(sequence):
        raise ValueError(
            f"Position {position} is out of range for sequence of "
            f"length {len(sequence)}."
        )

    if len(sequence) < 3:
        return ""

    if frame > 0:
        offset = frame - 1
        adjusted = position - 1 - offset
        if adjusted < 0:
            return ""
        codon_index = adjusted // 3
        codon_start = codon_index * 3 + offset
        codon = sequence[codon_start:codon_start + 3]
    else:
        rev_seq = reverse_complement(sequence)
        rev_pos = len(sequence) - position + 1
        offset = abs(frame) - 1
        adjusted = rev_pos - 1 - offset
        if adjusted < 0:
            return ""
        codon_index = adjusted // 3
        codon_start = codon_index * 3 + offset
        codon = rev_seq[codon_start:codon_start + 3]

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


def is_in_cds(position: int, cds_regions: list[tuple[int, int]]) -> bool:
    """Returns True if position falls within any CDS region.

    Args:
        position: 1-indexed SNP position.
        cds_regions: List of (start, end) tuples, both 1-indexed inclusive.

    Returns:
        bool: True if position is within at least one region.
    """
    return any(start <= position <= end for start, end in cds_regions)


def annotate_snp_with_regions(
    position: int,
    ref_sequence: str,
    alt_sequence: str,
    cds_regions: list[tuple[int, int]] | None = None,
    frame: int = 1,
) -> str:
    """Annotates a SNP respecting CDS boundaries and reading frame.

    If cds_regions is None, falls back to annotate_snp() (entire sequence
    treated as coding — backwards-compatible behavior).
    If cds_regions is an empty list, every position returns NON_CODING.
    Otherwise only positions within a CDS region receive functional
    annotation; positions outside return NON_CODING.

    Args:
        position: 1-indexed SNP position (forward strand).
        ref_sequence: Full reference DNA sequence (uppercase).
        alt_sequence: Full alternate DNA sequence (uppercase).
        cds_regions: List of (start, end) tuples, or None.
        frame: Reading frame. One of {1, 2, 3, -1, -2, -3}. Default 1.

    Returns:
        str: "SYNONYMOUS" | "NON_SYNONYMOUS" | "NONSENSE" | "NON_CODING"
    """
    if cds_regions is None:
        return annotate_snp(position, ref_sequence, alt_sequence, frame)

    if not is_in_cds(position, cds_regions):
        return "NON_CODING"

    return annotate_snp(position, ref_sequence, alt_sequence, frame)



def annotate_snp(
    position: int,
    ref_sequence: str,
    alt_sequence: str,
    frame: int = 1,
) -> str:
    """Annotates the functional effect of a SNP on the encoded amino acid.

    Identifies the codon containing the SNP in both reference and alternate
    sequences under the given reading frame, translates both, and classifies
    the change.

    Args:
        position: 1-indexed position of the SNP (forward strand).
        ref_sequence: Full reference DNA sequence (uppercase).
        alt_sequence: Full alternate (sample) DNA sequence (uppercase).
        frame: Reading frame. One of {1, 2, 3, -1, -2, -3}. Default 1.

    Returns:
        str: One of "SYNONYMOUS", "NON_SYNONYMOUS", "NONSENSE", "NON_CODING".

    Raises:
        ValueError: If frame is not one of {1, 2, 3, -1, -2, -3}.

    Note:
        Returns "NON_CODING" (instead of raising) when a sequence contains
        non-ACGT characters or an out-of-range position is encountered.
        Invalid frame values always raise ValueError — they represent a
        programming error, not a data quality issue.
    """
    _VALID_FRAMES = {1, 2, 3, -1, -2, -3}
    if frame not in _VALID_FRAMES:
        raise ValueError(
            f"Invalid frame '{frame}'. Must be one of {_VALID_FRAMES}."
        )

    try:
        ref_codon = get_codon(ref_sequence, position, frame)
        alt_codon = get_codon(alt_sequence, position, frame)
    except ValueError:
        return "NON_CODING"

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
