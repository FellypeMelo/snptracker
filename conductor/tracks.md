# Project Tracks

This file tracks all major tracks for the project. Each track has its own detailed plan in its respective folder.

---

## fasta_cli_20260227 — FASTA File Support & Argparse CLI
Folder: `conductor/tracks/fasta_cli_20260227/`
Status: ✅ Complete

## cds_regions_20260304 — CDS Regions (Coding vs Non-Coding)
Folder: N/A (single-session feature, no dedicated track folder)
Status: ✅ Complete

Changes:
- Added `is_in_cds()` and `annotate_snp_with_regions()` to `annotation.py`
- Added `parse_cds_regions()` to `main.py`
- `detect_snps()` gains optional `cds_regions` parameter (backwards-compatible)
- `parse_args()` gains `--cds` argument (`'start-end,start-end'` format)
- `_run_single_sample_mode` wires `--cds` into `detect_snps`
- 20 new tests added; 90 total passing, 94% coverage


Folder: N/A (single-session feature, no dedicated track folder)
Status: ✅ Complete

Changes:
- Added `get_trinucleotide_context()` to `main.py`
- Added `"context"` key to SNP dicts (SNPs only, not INDELs)
- Updated `print_snp_report` and `generate_snp_file` with Contexto column
- 8 new tests added; 70 total passing, 93% coverage
- Commit: `736f2be`

