# SNPTracker Technology Stack

## Core
- **Language:** Python 3.7+ (Standard library focus)
- **Frameworks:** Pure Python (No external dependencies)

## Modules
- **`main.py`:** CLI entrypoint, SNP detection, mutation classification, reporting.
- **`fasta_parser.py`:** FASTA file parsing (single and multi-sequence).
- **`annotation.py`:** Codon translation and SNP functional annotation
  (SYNONYMOUS / NON_SYNONYMOUS / NONSENSE / NON_CODING).
  Uses the complete standard genetic code (64 codons).
  Assumes reading frame starts at position 1 (frame +1).

## Data
- **Format:** FASTA files, plain text sequence inputs, text-based SNP reports.
- **Planned Support:** FASTQ input, VCF output.

## Future Infrastructure
- **Web Dashboard:** (To be determined, likely React or Streamlit).
- **External Integration:** Biopython (optional for complex parsing), requests (for database API access).
