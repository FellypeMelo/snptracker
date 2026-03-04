# Project Tracks

This file tracks all major tracks for the project. Each track has its own detailed plan in its respective folder.

---

## fasta_cli_20260227 — FASTA File Support & Argparse CLI
Folder: `conductor/tracks/fasta_cli_20260227/`
Status: ✅ Complete

## trinucleotide_context_20260304 — Trinucleotide Context (COSMIC format)
Folder: N/A (single-session feature, no dedicated track folder)
Status: ✅ Complete

Changes:
- Added `get_trinucleotide_context()` to `main.py`
- Added `"context"` key to SNP dicts (SNPs only, not INDELs)
- Updated `print_snp_report` and `generate_snp_file` with Contexto column
- 8 new tests added; 70 total passing, 93% coverage
- Commit: `736f2be`

