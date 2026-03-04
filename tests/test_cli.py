import unittest
import os
from main import parse_args, load_sequence, run_multi_sample
from fasta_parser import read_all_sequences

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


class TestMultiSampleCLI(unittest.TestCase):

    def setUp(self):
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            if os.path.exists(f):
                os.remove(f)

    def create_temp_fasta(self, content):
        file_path = f"cli_multi_temp_{len(self.temp_files)}.fasta"
        with open(file_path, 'w') as f:
            f.write(content)
        self.temp_files.append(file_path)
        return file_path

    def test_parse_args_input_flag(self):
        """Test that --input argument is accepted."""
        args = parse_args(["--input", "data/sequences.txt"])
        self.assertEqual(args.input, "data/sequences.txt")

    def test_parse_args_input_mutually_exclusive_with_reference(self):
        """Test that --input and --reference cannot be used together."""
        with self.assertRaises(SystemExit):
            parse_args(["--input", "file.fasta", "--reference", "ACTG", "--sample", "ACTT"])

    def test_parse_args_input_no_sample_required(self):
        """Test that --sample is not required when --input is used."""
        args = parse_args(["--input", "data/sequences.txt"])
        self.assertIsNone(args.reference)
        self.assertIsNone(args.sample)

    def test_multi_sample_integration(self):
        """Test full flow: read multi-sequence FASTA and run analysis."""
        path = self.create_temp_fasta(
            ">reference\nACTG\n>sample1\nACTT\n>sample2\nACTG"
        )
        sequences = read_all_sequences(path)
        ref_header, ref_seq = sequences[0]
        samples = sequences[1:]
        result = run_multi_sample(ref_seq, samples)
        self.assertIn("sample1", result)
        self.assertIn("sample2", result)
        self.assertEqual(len(result["sample1"]), 1)
        self.assertEqual(len(result["sample2"]), 0)


    def test_parse_args_cds_single_region(self):
        """--cds '1-90' is accepted and stored."""
        args = parse_args(["--reference", "ACTG", "--sample", "ACTT", "--cds", "1-90"])
        self.assertEqual(args.cds, "1-90")

    def test_parse_args_cds_multiple_regions(self):
        """--cds '1-90,100-150' is accepted."""
        args = parse_args(["--reference", "ACTG", "--sample", "ACTT", "--cds", "1-90,100-150"])
        self.assertEqual(args.cds, "1-90,100-150")

    def test_parse_args_cds_default_is_none(self):
        """--cds is optional; defaults to None."""
        args = parse_args(["--reference", "ACTG", "--sample", "ACTT"])
        self.assertIsNone(args.cds)


if __name__ == "__main__":
    unittest.main()
