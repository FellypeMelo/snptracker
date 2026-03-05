"""
SNPTracker - Functional Prediction Module

Provides Grantham Score-based functional impact prediction for
NON_SYNONYMOUS SNPs using physicochemical amino acid properties
(composition, polarity, volume) as described in:

    Grantham R. (1974). Amino acid difference formula to help
    explain protein evolution. Science 185(4154):862-864.

Classification thresholds:
    CONSERVATIVE  — score  0–50   (chemically similar substitution)
    MODERATE      — score 51–100  (moderately different substitution)
    RADICAL       — score  >100   (chemically dissimilar substitution)

Future integration intent:
    SIFT (sift.bii.a-star.edu.sg) and PolyPhen-2 external APIs are
    planned for a future milestone. Current blockers:
    - SIFT direct API uses asynchronous batch processing (10–15 min),
      impractical for CLI integration.
    - PolyPhen-2 public server is currently unavailable.
    - Ensembl VEP REST API requires chromosomal coordinates and a
      species/gene identifier not available from raw sequences alone.
    When the above blockers are resolved, query_sift() and
    query_polyphen() can be added here using urllib (stdlib only).

Public API:
    translate_dna_to_protein(sequence, frame=1) -> list[str]
    get_amino_acid_change(snp, ref_sequence, frame=1)
        -> tuple[str, int, str] | None
    grantham_score(ref_aa, alt_aa) -> int | None
    grantham_prediction(score) -> str
    predict_functional_impact(snp, ref_sequence, frame=1) -> dict
"""

import math
from annotation import get_codon, translate_codon, reverse_complement

# ---------------------------------------------------------------------------
# Grantham (1974) physicochemical properties: (composition, polarity, volume)
# ---------------------------------------------------------------------------
_GRANTHAM_PROPERTIES: dict[str, tuple[float, float, float]] = {
    "Ala": (0,    8.1,  31.0),
    "Arg": (0.65, 10.5, 124.0),
    "Asn": (1.33, 11.6,  56.0),
    "Asp": (1.38, 13.0,  54.0),
    "Cys": (2.75,  5.5,  55.0),
    "Gln": (0.89, 10.5,  85.0),
    "Glu": (0.92, 12.3,  83.0),
    "Gly": (0.74,  9.0,   3.0),
    "His": (0.58, 10.4,  96.0),
    "Ile": (0,     5.2, 111.0),
    "Leu": (0,     4.9, 111.0),
    "Lys": (0.33, 11.3, 119.0),
    "Met": (0,     5.7, 105.0),
    "Phe": (0,     5.2, 132.0),
    "Pro": (0.39,  8.0,  32.5),
    "Ser": (1.42,  9.2,  32.0),
    "Thr": (0.71,  8.6,  61.0),
    "Trp": (0.13,  5.4, 170.0),
    "Tyr": (0.20,  6.2, 136.0),
    "Val": (0,     5.9,  84.0),
}

# Normalization coefficients (empirically derived from known reference pairs)
_C_COEF = 4736.0
_P_COEF = 358.8
_V_COEF = 0.912


def translate_dna_to_protein(sequence: str, frame: int = 1) -> list[str]:
    """Translates a DNA sequence to a list of amino acid 3-letter codes.

    Uses the existing CODON_TABLE from annotation.py. Translation stops
    at the first STOP codon (inclusive, represented as '*'). Incomplete
    trailing codons are ignored.

    Args:
        sequence: Full DNA sequence string (uppercase).
        frame: Reading frame. One of {1, 2, 3, -1, -2, -3}. Default 1.

    Returns:
        list[str]: Amino acid codes (e.g. ["Met", "Val"]) with "*" for
            STOP codons. Empty list if the sequence yields no complete
            codons in the given frame.
    """
    if frame < 0:
        seq = reverse_complement(sequence)
        offset = abs(frame) - 1
    else:
        seq = sequence
        offset = frame - 1

    protein = []
    i = offset
    while i + 3 <= len(seq):
        codon = seq[i:i + 3]
        aa = translate_codon(codon)
        if aa == "STOP":
            protein.append("*")
            break
        protein.append(aa)
        i += 3

    return protein


def get_amino_acid_change(
    snp: dict,
    ref_sequence: str,
    frame: int = 1,
) -> tuple[str, int, str] | None:
    """Returns the amino acid change for a NON_SYNONYMOUS SNP.

    Constructs the alternate sequence from the SNP dict, retrieves both
    codons via get_codon(), and returns the amino acids together with
    the 1-indexed amino acid position in the protein.

    Args:
        snp: SNP dict as produced by detect_snps(). Must contain keys
            'position', 'alternate', and 'annotation'.
        ref_sequence: Full reference DNA sequence (uppercase).
        frame: Reading frame. One of {1, 2, 3, -1, -2, -3}. Default 1.

    Returns:
        tuple[str, int, str]: (ref_aa, aa_position, alt_aa) when
            snp['annotation'] == 'NON_SYNONYMOUS'.
        None: for SYNONYMOUS, NONSENSE, NON_CODING, or INDEL variants.
    """
    if snp.get("annotation") != "NON_SYNONYMOUS":
        return None

    position = snp["position"]
    alt_base = snp["alternate"]

    alt_sequence = (
        ref_sequence[:position - 1] + alt_base + ref_sequence[position:]
    )

    ref_codon = get_codon(ref_sequence, position, frame)
    if not ref_codon:
        return None

    ref_aa = translate_codon(ref_codon)

    if frame > 0:
        offset = frame - 1
        aa_position = (position - 1 - offset) // 3 + 1
    else:
        rev_pos = len(ref_sequence) - position + 1
        offset = abs(frame) - 1
        aa_position = (rev_pos - 1 - offset) // 3 + 1

    alt_codon = get_codon(alt_sequence, position, frame)
    alt_aa = translate_codon(alt_codon)

    return (ref_aa, aa_position, alt_aa)


def grantham_score(ref_aa: str, alt_aa: str) -> int | None:
    """Computes the Grantham distance between two amino acids.

    Uses physicochemical properties (composition, polarity, volume)
    from Grantham (1974). Returns None if either amino acid is not a
    standard residue (e.g., "STOP" or unknown).

    Args:
        ref_aa: Reference amino acid (3-letter code, e.g. "Val").
        alt_aa: Alternate amino acid (3-letter code, e.g. "Ile").

    Returns:
        int: Grantham distance (0 = identical, higher = more dissimilar).
        None: If either amino acid is not in the property table.
    """
    if ref_aa not in _GRANTHAM_PROPERTIES or alt_aa not in _GRANTHAM_PROPERTIES:
        return None

    c1, p1, v1 = _GRANTHAM_PROPERTIES[ref_aa]
    c2, p2, v2 = _GRANTHAM_PROPERTIES[alt_aa]

    dc, dp, dv = c1 - c2, p1 - p2, v1 - v2
    return round(math.sqrt(
        _C_COEF * dc ** 2 + _P_COEF * dp ** 2 + _V_COEF * dv ** 2
    ))


def grantham_prediction(score: int) -> str:
    """Classifies a Grantham score into a qualitative category.

    Args:
        score: Non-negative Grantham distance value.

    Returns:
        str: "CONSERVATIVE" (0–50), "MODERATE" (51–100), or
             "RADICAL" (>100).
    """
    if score <= 50:
        return "CONSERVATIVE"
    if score <= 100:
        return "MODERATE"
    return "RADICAL"


def predict_functional_impact(
    snp: dict,
    ref_sequence: str,
    frame: int = 1,
) -> dict:
    """Predicts the functional impact of a SNP using Grantham Score.

    Only applicable to NON_SYNONYMOUS SNPs. All other annotation
    categories (SYNONYMOUS, NONSENSE, NON_CODING) and INDELs return
    an empty dict.

    Args:
        snp: SNP dict as produced by detect_snps(). Must contain
            'position', 'alternate', and 'annotation'.
        ref_sequence: Full reference DNA sequence (uppercase).
        frame: Reading frame. One of {1, 2, 3, -1, -2, -3}. Default 1.

    Returns:
        dict: {'grantham_score': int, 'grantham_prediction': str}
              when the SNP is NON_SYNONYMOUS and amino acids are known.
              {} for any other variant type.
    """
    change = get_amino_acid_change(snp, ref_sequence, frame)
    if change is None:
        return {}

    ref_aa, _aa_pos, alt_aa = change
    score = grantham_score(ref_aa, alt_aa)
    if score is None:
        return {}

    return {
        "grantham_score": score,
        "grantham_prediction": grantham_prediction(score),
    }
