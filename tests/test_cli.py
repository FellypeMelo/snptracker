import unittest
from main import parse_args

class TestCLI(unittest.TestCase):
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
