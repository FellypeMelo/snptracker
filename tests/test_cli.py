import unittest
import os
from main import parse_args, load_sequence

class TestCLI(unittest.TestCase):
    def setUp(self):
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            if os.path.exists(f):
                os.remove(f)

    def create_temp_fasta(self, content):
        file_path = f"cli_temp_{len(self.temp_files)}.fasta"
        with open(file_path, 'w') as f:
            f.write(content)
        self.temp_files.append(file_path)
        return file_path

    def test_load_sequence_from_string(self):
        """Test that a raw string is returned as-is if no file exists."""
        self.assertEqual(load_sequence("ACTG"), "ACTG")

    def test_load_sequence_from_file(self):
        """Test that a sequence is loaded from a FASTA file."""
        path = self.create_temp_fasta(">seq\nACTG")
        self.assertEqual(load_sequence(path), "ACTG")

    def test_parse_args_valid(self):
        """Test that valid arguments are parsed correctly."""
        args = ["--reference", "ACTG", "--sample", "ACTT", "--output", "test_report.txt"]
        parsed = parse_args(args)
        self.assertEqual(parsed.reference, "ACTG")
        self.assertEqual(parsed.sample, "ACTT")
        self.assertEqual(parsed.output, "test_report.txt")

    def test_parse_args_default_output(self):
        """Test that the default output file is used if not provided."""
        args = ["--reference", "ACTG", "--sample", "ACTT"]
        parsed = parse_args(args)
        self.assertEqual(parsed.output, "snps_report.txt")

    def test_parse_args_missing_required(self):
        """Test that missing required arguments raise an error (or exit)."""
        # argparse usually exits with status 2 on missing required args
        with self.assertRaises(SystemExit):
            parse_args(["--reference", "ACTG"])

if __name__ == "__main__":
    unittest.main()
