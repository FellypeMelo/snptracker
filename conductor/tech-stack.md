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

- **`fasta_parser.py`:** FASTA file parsing (single and multi-sequence).
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
}
```

## Data
- **Format:** FASTA files, plain text sequence inputs, text-based SNP reports.
- **Planned Support:** FASTQ input, VCF output.

## Future Infrastructure
- **Web Dashboard:** (To be determined, likely React or Streamlit).
- **External Integration:** Biopython (optional for complex parsing), requests (for database API access).
