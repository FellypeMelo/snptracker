# Task Notes: Set up CLI argument parsing with argparse
... (existing notes) ...

# Task Notes: Implement FASTA file parser
... (existing notes) ...

# Task Notes: Integrate FASTA parser into CLI

## Scope
- **Goal:** Allow users to provide file paths instead of raw sequence strings in the CLI.
- **Files to Change:**
    - `main.py`: Update to detect if inputs are files and use the parser.
    - `tests/test_cli.py`: Add integration tests for file-based inputs.
- **Expected Behavior:**
    - `python main.py --reference ref.fasta --sample smp.fasta` should work.
    - If the argument is a valid file path ending in `.fasta` or `.fa`, it should be parsed.
    - If it's not a file path, it should still support raw strings (for backward compatibility).
- **Edge Cases:**
    - File exists but is not a valid FASTA.
    - File does not exist (should probably error out or treat as raw string - I'll implement it to check if it's a file first).
- **What should NOT change:**
    - Core SNP detection logic.

## Proposed Implementation
- **Logic in `main.py`:**
    1. Define a helper `get_sequence(input_str)`:
        - If `os.path.isfile(input_str)`: return `read_fasta(input_str)`.
        - Else: return `input_str`.
    2. In `main()`, call `get_sequence` for both reference and sample.
- **Testing Strategy:**
    - Create temporary FASTA files in tests.
    - Run `parse_args` and then verify the resolved sequences.
