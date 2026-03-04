# Specification: Multi-Sample Support via Single FASTA File

## Goal
Extend SNPTracker to detect SNPs across multiple samples by reading a single
multi-sequence FASTA file, where the first sequence is the reference and all
remaining sequences are samples.

## Scope
- **FASTA Parser:** New `read_all_sequences()` function returning all sequences.
- **Logic:** New `run_multi_sample()` function orchestrating detection per sample.
- **CLI:** New `--input` flag (mutually exclusive with `--reference`).

## Requirements
- First sequence in the file is always the reference.
- All subsequent sequences are treated as samples.
- One report file generated per sample with SNPs (samples with 0 SNPs skip file).
- `--reference`/`--sample` mode continues to work unchanged.
