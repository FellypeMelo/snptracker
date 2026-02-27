# Task Notes: Set up CLI argument parsing with argparse

## Scope
- **Goal:** Replace hardcoded sequences in `main.py` with CLI arguments using `argparse`.
- **Files to Change:**
    - `main.py`: Refactor to use `argparse`.
- **New Files:**
    - `tests/test_cli.py`: Unit tests for CLI parsing.
- **Expected Behavior:**
    - Running `python main.py --reference <seq1> --sample <seq2>` should output the SNP report.
    - Providing invalid arguments should show a help message.
- **Edge Cases:**
    - Missing required arguments.
    - Empty sequences.
- **What should NOT change:**
    - Core SNP detection logic in `detect_snps` and `classify_mutation`.
    - Report formatting logic in `print_snp_report`.

## Proposed Implementation
- **Changes in `main.py`:**
    - Add `import argparse`.
    - Create `parse_args(args=None)` function.
    - Update `main()` to use `parse_args()`.
- **Testing Strategy:**
    - Use `unittest` or `pytest` to verify `parse_args` with various input lists.
    - Mock `sys.stdout` if needed to verify report output in later tasks.
