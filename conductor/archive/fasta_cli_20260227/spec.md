# Specification: Implement FASTA file support and Argparse CLI

## Goal
Enhance SNPTracker to support standard FASTA file inputs and provide a robust Command Line Interface (CLI) using Python's `argparse`.

## Scope
- **CLI Enhancement:** Replace manual variable editing in `main.py` with an `argparse` interface supporting `--reference` and `--sample` arguments.
- **FASTA Parsing:** Implement a parser to read DNA sequences from FASTA files.
- **Integration:** Update the core detection logic to handle sequences loaded from files.

## Requirements
- Support multi-line FASTA files.
- Handle basic FASTA headers (starting with '>').
- Provide helpful CLI usage messages and error handling for missing files.
