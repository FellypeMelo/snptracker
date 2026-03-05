# SNPTracker Technology Stack

## Core
- **Language:** Python 3.7+ (Standard library focus)
- **Frameworks:** Pure Python (No external dependencies)

## Modules
- **`main.py`:** CLI entrypoint, SNP detection, mutation classification, reporting.
  Key functions:
  - `detect_snps(reference, sample, cds_regions=None, frame=1)` — compares two sequences base-by-base and returns a list of SNP dicts. Accepts optional `cds_regions` to restrict functional annotation to coding regions, and `frame` to select the reading frame (1/2/3/-1/-2/-3).
  - `classify_mutation(ref_base, alt_base)` — classifies TRANSITION or TRANSVERSION.
  - `parse_cds_regions(cds_str)` — parses CLI string `'1-90,100-150'` into `[(1,90),(100,150)]`.
  - `get_trinucleotide_context(reference, position, ref_base, alt_base)` — returns COSMIC-format trinucleotide context (`X[R>A]Y`).
  - `run_multi_sample(reference, samples)` — batch detection across multiple samples.
  - `print_snp_report(snps, reference, sample, frame=1)` / `generate_snp_file` — output formatting (terminal and file). Report header includes the active reading frame.
  - `load_sequence(input_data)` — returns a raw sequence string or reads from a FASTA file. Raises `FileNotFoundError` if the argument looks like a file path (has an extension or path separator) but the file does not exist, preventing silent mis-annotation from typos.
  - `parse_args` / `main` — CLI wiring. Supports `--frame` argument (choices: 1, 2, 3, -1, -2, -3; default: 1).

- **`prediction.py`:** Functional impact prediction for NON_SYNONYMOUS SNPs.
  Uses Grantham Score (1974) based on amino acid physicochemical properties
  (composition, polarity, volume). No external dependencies required.
  Key functions:
  - `translate_dna_to_protein(sequence, frame=1)` — translates a DNA sequence to a list of 3-letter amino acid codes. Supports all six reading frames; stops at STOP codon (included as `"*"`). Incomplete trailing codons are ignored.
  - `get_amino_acid_change(snp, ref_sequence, frame=1)` — returns `(ref_aa, aa_position, alt_aa)` for NON_SYNONYMOUS SNPs; `None` for all other types (SYNONYMOUS, NONSENSE, NON_CODING, INDELs).
  - `grantham_score(ref_aa, alt_aa)` — computes Grantham distance between two 3-letter amino acid codes. Returns `None` for unknown residues (e.g., `"STOP"`).
  - `grantham_prediction(score)` — classifies score: `CONSERVATIVE` (0–50), `MODERATE` (51–100), `RADICAL` (>100).
  - `predict_functional_impact(snp, ref_sequence, frame=1)` — entry point. Returns `{'grantham_score': int, 'grantham_prediction': str}` for NON_SYNONYMOUS SNPs; `{}` for all others.


  Key functions:
  - `read_fasta(file_path)` — returns the first sequence as a plain string.
  - `read_all_sequences(file_path)` — returns all sequences as `list[tuple[header, sequence]]`.

- **`annotation.py`:** Codon translation and SNP functional annotation.
  Uses the complete standard genetic code (64 codons).
  Supports six reading frames: +1, +2, +3 (forward strand) and -1, -2, -3
  (reverse complement). Default frame is +1.
  Key functions:
  - `reverse_complement(sequence)` — returns the reverse complement of a DNA string; raises `ValueError` for non-ACGT input.
  - `get_codon(sequence, position, frame=1)` — returns the triplet containing the given 1-indexed position under the specified frame. Returns `""` for positions before the frame start or in incomplete trailing codons.
  - `translate_codon(codon)` — returns amino acid abbreviation or `"STOP"`.
  - `annotate_snp(position, ref_sequence, alt_sequence, frame=1)` — returns one of the four annotation strings (entire sequence treated as coding). Raises `ValueError` for invalid frame values; returns `NON_CODING` for non-ACGT sequences or out-of-range positions.
  - `is_in_cds(position, cds_regions)` — returns True if position falls within any CDS region.
  - `annotate_snp_with_regions(position, ref_sequence, alt_sequence, cds_regions=None, frame=1)` — like `annotate_snp` but respects CDS boundaries; positions outside any region return `NON_CODING`.

## SNP Dict Format

Each detected SNP is represented as a Python dict. Fields:

```python
{
    "position":  int,   # 1-indexed position in the sequence
    "reference": str,   # Reference base (A/C/G/T or '-' for insertions)
    "alternate": str,   # Alternate base (A/C/G/T or '-' for deletions)
    "type":      str,   # "TRANSITION" | "TRANSVERSION" | "INSERTION" | "DELETION"
    "annotation":str,   # "SYNONYMOUS" | "NON_SYNONYMOUS" | "NONSENSE" | "NON_CODING"
    "context":   str,   # COSMIC trinucleotide format, e.g. "T[A>G]C" — SNPs only
                        # Not present in INSERTION / DELETION dicts
    # --- Optional: present only when --predict is active and annotation == NON_SYNONYMOUS ---
    "grantham_score":      int,  # Grantham distance (0 = identical, higher = more dissimilar)
    "grantham_prediction": str,  # "CONSERVATIVE" | "MODERATE" | "RADICAL"
}
```

## Data
- **Format:** FASTA files, plain text sequence inputs, text-based SNP reports.
- **Planned Support:** FASTQ input, VCF output.

## Future Infrastructure
- **Web Dashboard:** (To be determined, likely React or Streamlit).
- **External Integration:** Biopython (optional for complex parsing), requests (for database API access).

## Functional Prediction — Future API Integration (Planned)

The current implementation uses **Grantham Score** (local, no network).
Integration with external prediction tools is planned for a future milestone.

**Intended tools:** SIFT and PolyPhen-2.

**Current blockers:**
- SIFT (`sift.bii.a-star.edu.sg`): Uses web-form-based asynchronous batch
  processing (~10–15 min per submission). Not suitable for synchronous CLI use.
  A polling-with-retry pattern would be required.
- PolyPhen-2 (`genetics.bwh.harvard.edu/pph2`): Public server currently
  returns 404. Needs a stable alternative endpoint.
- Ensembl VEP REST API: Provides SIFT + PolyPhen scores but requires
  chromosomal coordinates and a species/gene identifier — information not
  available from raw DNA sequences alone.

**Implementation path when blockers are resolved:**
1. Add `query_sift(protein_sequence, position, ref_aa, alt_aa, timeout)` to `prediction.py` using `urllib` (stdlib).
2. Add `query_polyphen(...)` similarly.
3. Extend `predict_functional_impact()` to call both APIs and merge results into the SNP dict alongside the existing Grantham score.
4. Grantham Score remains as a fast local fallback when APIs are unavailable.

