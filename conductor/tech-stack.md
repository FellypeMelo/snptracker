# SNPTracker Technology Stack

## Core
- **Language:** Python 3.7+ (Standard library focus)
- **Frameworks:** Pure Python (No external dependencies)

## Modules
- **`main.py`:** CLI entrypoint, SNP detection, mutation classification, reporting.
  Key functions:
  - `detect_snps(reference, sample)` — compares two sequences base-by-base and returns a list of SNP dicts.
  - `classify_mutation(ref_base, alt_base)` — classifies TRANSITION or TRANSVERSION.
  - `get_trinucleotide_context(reference, position, ref_base, alt_base)` — returns COSMIC-format trinucleotide context (`X[R>A]Y`).
  - `run_multi_sample(reference, samples)` — batch detection across multiple samples.
  - `print_snp_report` / `generate_snp_file` — output formatting (terminal and file).
  - `load_sequence` / `parse_args` / `main` — CLI wiring.

- **`fasta_parser.py`:** FASTA file parsing (single and multi-sequence).
  Key functions:
  - `read_fasta(file_path)` — returns the first sequence as a plain string.
  - `read_all_sequences(file_path)` — returns all sequences as `list[tuple[header, sequence]]`.

- **`annotation.py`:** Codon translation and SNP functional annotation.
  Uses the complete standard genetic code (64 codons).
  Assumes reading frame starts at position 1 (frame +1).
  Key functions:
  - `get_codon(sequence, position)` — returns the triplet containing the given 1-indexed position.
  - `translate_codon(codon)` — returns amino acid abbreviation or `"STOP"`.
  - `annotate_snp(position, ref_sequence, alt_sequence)` — returns `"SYNONYMOUS"`, `"NON_SYNONYMOUS"`, `"NONSENSE"`, or `"NON_CODING"`.

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
