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


## reading_frames_20260305 — Multiple Reading Frames
Folder: N/A (single-session feature, no dedicated track folder)
Status: ✅ Complete

Changes:
- Added `reverse_complement()` to `annotation.py`
- `get_codon()` gains `frame` param: +1/+2/+3 (forward), -1/-2/-3 (reverse complement)
- `annotate_snp()` and `annotate_snp_with_regions()` gain `frame` param (default=1)
- `detect_snps()` gains `frame` param (default=1, backwards-compatible)
- `parse_args()` gains `--frame` argument (choices: 1,2,3,-1,-2,-3; default: 1)
- `print_snp_report()` displays active reading frame in header
- Commit: `52b46e1`

## bugfix_frame_robustness_20260305 — Frame Robustness & File Path Validation
Folder: N/A
Status: ✅ Complete

Changes:
- Fix 1 (annotation.py): `annotate_snp()` validates `frame` before codon lookup,
  then wraps `get_codon()` in try/except ValueError. Invalid-base sequences
  (e.g., a file path passed as a sequence) return NON_CODING instead of crashing.
  Invalid frame values still raise ValueError.
- Fix 2 (main.py): `load_sequence()` raises FileNotFoundError when the argument
  has a file extension or directory separator but the file does not exist.
  Prevents a mistyped filename from being silently treated as a raw DNA sequence.
- Root cause: reverse frames call reverse_complement() which validates characters
  strictly; forward frames swallowed the same error silently via translate_codon
  try/except. The asymmetry is now resolved.
- 13 new regression tests; 151 total passing, 95% coverage.
