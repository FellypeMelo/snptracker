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


class TestParseArgsFrame(unittest.TestCase):

    def test_frame_default_is_1(self):
        """--frame is optional; defaults to 1."""
        args = parse_args(["--reference", "ACTG", "--sample", "ACTT"])
        self.assertEqual(args.frame, 1)

    def test_frame_2_accepted(self):
        """--frame 2 is accepted and stored as int."""
        args = parse_args(
            ["--reference", "ACTG", "--sample", "ACTT", "--frame", "2"]
        )
        self.assertEqual(args.frame, 2)

    def test_frame_3_accepted(self):
        args = parse_args(
            ["--reference", "ACTG", "--sample", "ACTT", "--frame", "3"]
        )
        self.assertEqual(args.frame, 3)

    def test_frame_minus1_accepted(self):
        """--frame -1 is accepted."""
        args = parse_args(
            ["--reference", "ACTG", "--sample", "ACTT", "--frame", "-1"]
        )
        self.assertEqual(args.frame, -1)

    def test_frame_minus3_accepted(self):
        args = parse_args(
            ["--reference", "ACTG", "--sample", "ACTT", "--frame", "-3"]
        )
        self.assertEqual(args.frame, -3)

    def test_frame_invalid_value_raises(self):
        """--frame 0 is not a valid frame and causes SystemExit."""
        with self.assertRaises(SystemExit):
            parse_args(
                ["--reference", "ACTG", "--sample", "ACTT", "--frame", "0"]
            )

    def test_frame_invalid_string_raises(self):
        """--frame abc causes SystemExit (not an int)."""
        with self.assertRaises(SystemExit):
            parse_args(
                ["--reference", "ACTG", "--sample", "ACTT", "--frame", "abc"]
            )


class TestReportWithFrame(unittest.TestCase):

    def test_report_header_shows_frame(self):
        """print_snp_report includes the active frame in the header."""
        from unittest.mock import patch
        import io
        from main import print_snp_report
        snps = []
        with patch("sys.stdout", new_callable=io.StringIO) as mock_out:
            print_snp_report(snps, "ACTG", "ACTG", frame=2)
            output = mock_out.getvalue()
        self.assertIn("Frame", output)
        self.assertIn("+2", output)

    def test_report_header_frame1_shows_plus1(self):
        """Default frame=1 is displayed as +1."""
        from unittest.mock import patch
        import io
        from main import print_snp_report
        with patch("sys.stdout", new_callable=io.StringIO) as mock_out:
            print_snp_report([], "ACTG", "ACTG", frame=1)
            output = mock_out.getvalue()
        self.assertIn("+1", output)

    def test_report_header_negative_frame(self):
        """Negative frame is displayed with minus sign."""
        from unittest.mock import patch
        import io
        from main import print_snp_report
        with patch("sys.stdout", new_callable=io.StringIO) as mock_out:
            print_snp_report([], "ACTG", "ACTG", frame=-1)
            output = mock_out.getvalue()
        self.assertIn("-1", output)


class TestLoadSequencePathValidation(unittest.TestCase):
    """Tests for load_sequence() file-path detection and error reporting."""

    def test_raw_acgt_string_returned_as_is(self):
        """A plain ACGT sequence string is returned unchanged."""
        self.assertEqual(load_sequence("ACTG"), "ACTG")

    def test_existing_fasta_file_loads_sequence(self):
        """An existing .fasta file is parsed and its sequence returned."""
        import tempfile, os
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(">seq\nACTG")
            path = f.name
        try:
            self.assertEqual(load_sequence(path), "ACTG")
        finally:
            os.unlink(path)

    def test_missing_fasta_extension_raises_file_not_found(self):
        """A non-existent path with .fasta extension raises FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            load_sequence("nonexistent_file.fasta")

    def test_missing_fa_extension_raises_file_not_found(self):
        """A non-existent path with .fa extension raises FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            load_sequence("ref.fa")

    def test_missing_txt_extension_raises_file_not_found(self):
        """A non-existent path with .txt extension raises FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            load_sequence("sequences.txt")

    def test_path_with_directory_separator_raises_file_not_found(self):
        """A path containing a directory separator raises FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            load_sequence("data/ref.fasta")

    def test_absolute_missing_path_raises_file_not_found(self):
        """An absolute path that does not exist raises FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            load_sequence("/tmp/does_not_exist_snptracker.fasta")

    def test_error_message_contains_path(self):
        """FileNotFoundError message includes the offending path."""
        path = "missing_ref.fasta"
        with self.assertRaises(FileNotFoundError) as ctx:
            load_sequence(path)
        self.assertIn(path, str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
