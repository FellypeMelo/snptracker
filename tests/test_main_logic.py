import unittest
import os
from unittest.mock import patch
import io
from main import detect_snps, classify_mutation, generate_snp_file, print_snp_report, run_multi_sample

class TestMainLogic(unittest.TestCase):
    def setUp(self):
        self.test_output = "test_output.txt"

    def tearDown(self):
        if os.path.exists(self.test_output):
            os.remove(self.test_output)

    def test_print_snp_report(self):
        """Test that the SNP report is printed correctly."""
        snps = [
            {'position': 1, 'reference': 'A', 'alternate': 'G', 'type': 'TRANSITION'}
        ]
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            print_snp_report(snps, "ACTG", "GCTG")
            output = fake_out.getvalue()
            self.assertIn("SNPTracker - Relatório de Mutações", output)
            self.assertIn("Total de variações encontradas: 1", output)
            self.assertIn("TRANSITION", output)

    def test_generate_snp_file(self):
        """Test that the SNP report file is generated correctly."""
        snps = [
            {'position': 1, 'reference': 'A', 'alternate': 'G', 'type': 'TRANSITION'}
        ]
        generate_snp_file(snps, output_file=self.test_output)
        self.assertTrue(os.path.exists(self.test_output))
        with open(self.test_output, 'r') as f:
            content = f.read()
            self.assertIn("SNPTRACKER - RELATÓRIO DE SNPs", content)
            self.assertIn("TRANSITION", content)

    def test_detect_snps_identical(self):
        """Test detection with identical sequences."""
        ref = "ACTG"
        smp = "ACTG"
        snps = detect_snps(ref, smp)
        self.assertEqual(len(snps), 0)

    def test_detect_snps_single_snp(self):
        """Test detection with a single SNP."""
        ref = "ACTG"
        smp = "ACTT"
        snps = detect_snps(ref, smp)
        self.assertEqual(len(snps), 1)
        self.assertEqual(snps[0]['position'], 4)
        self.assertEqual(snps[0]['reference'], 'G')
        self.assertEqual(snps[0]['alternate'], 'T')
        self.assertEqual(snps[0]['type'], 'TRANSVERSION')

    def test_detect_snps_includes_annotation_field(self):
        """Every SNP dict must contain the 'annotation' key."""
        snps = detect_snps("GCTAAA", "GTTAAA")
        self.assertTrue(len(snps) > 0)
        for snp in snps:
            self.assertIn('annotation', snp)

    def test_detect_snps_indel_annotation_is_non_coding(self):
        """INSERTION and DELETION SNPs must have annotation='NON_CODING'."""
        ins_snps = detect_snps("ACTG", "ACTGA")
        self.assertEqual(ins_snps[0]['annotation'], 'NON_CODING')

        del_snps = detect_snps("ACTG", "ACT")
        self.assertEqual(del_snps[0]['annotation'], 'NON_CODING')

    def test_detect_snps_insertion(self):
        """Test detection with an insertion."""
        ref = "ACTG"
        smp = "ACTGA"
        snps = detect_snps(ref, smp)
        self.assertEqual(len(snps), 1)
        self.assertEqual(snps[0]['type'], 'INSERTION')
        self.assertEqual(snps[0]['alternate'], 'A')

    def test_detect_snps_deletion(self):
        """Test detection with a deletion."""
        ref = "ACTG"
        smp = "ACT"
        snps = detect_snps(ref, smp)
        self.assertEqual(len(snps), 1)
        self.assertEqual(snps[0]['type'], 'DELETION')
        self.assertEqual(snps[0]['reference'], 'G')

    def test_classify_mutation_transition(self):
        """Test classification of transitions."""
        self.assertEqual(classify_mutation('A', 'G'), 'TRANSITION')
        self.assertEqual(classify_mutation('C', 'T'), 'TRANSITION')

    def test_classify_mutation_transversion(self):
        """Test classification of transversions."""
        self.assertEqual(classify_mutation('A', 'C'), 'TRANSVERSION')
        self.assertEqual(classify_mutation('G', 'T'), 'TRANSVERSION')


class TestRunMultiSample(unittest.TestCase):

    def test_run_multi_sample_empty_list(self):
        """Test that empty samples list returns empty dict."""
        result = run_multi_sample("ACTG", [])
        self.assertEqual(result, {})

    def test_run_multi_sample_single(self):
        """Test with one sample that has one SNP."""
        result = run_multi_sample("ACTG", [("sample1", "ACTT")])
        self.assertIn("sample1", result)
        self.assertEqual(len(result["sample1"]), 1)
        self.assertEqual(result["sample1"][0]["position"], 4)

    def test_run_multi_sample_multiple(self):
        """Test with N samples returns results for each."""
        samples = [("s1", "ACTT"), ("s2", "GCTG"), ("s3", "ACTG")]
        result = run_multi_sample("ACTG", samples)
        self.assertEqual(len(result), 3)
        self.assertEqual(len(result["s1"]), 1)   # 1 SNP
        self.assertEqual(len(result["s2"]), 1)   # 1 SNP
        self.assertEqual(len(result["s3"]), 0)   # idêntica

    def test_run_multi_sample_identical(self):
        """Test that sample identical to reference returns 0 SNPs."""
        result = run_multi_sample("ACTG", [("ref_copy", "ACTG")])
        self.assertEqual(result["ref_copy"], [])

    def test_run_multi_sample_with_indel(self):
        """Test that indels are detected per sample."""
        result = run_multi_sample("ACTG", [("del_sample", "ACT")])
        self.assertEqual(result["del_sample"][0]["type"], "DELETION")


class TestRunModes(unittest.TestCase):

    def setUp(self):
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            if os.path.exists(f):
                os.remove(f)

    def create_temp_fasta(self, content):
        path = f"mode_temp_{len(self.temp_files)}.fasta"
        with open(path, 'w') as f:
            f.write(content)
        self.temp_files.append(path)
        return path

    def test_single_sample_mode_prints_report(self):
        """Test _run_single_sample_mode produces output."""
        from main import _run_single_sample_mode, parse_args
        args = parse_args(["--reference", "ACTG", "--sample", "ACTT"])
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            _run_single_sample_mode(args)
            output = fake_out.getvalue()
        self.assertIn("SNPTracker", output)

    def test_multi_sample_mode_no_sequences(self):
        """Test _run_multi_sample_mode with empty file shows error."""
        from main import _run_multi_sample_mode, parse_args
        path = self.create_temp_fasta("")
        self.temp_files.append(path)
        args = parse_args(["--input", path])
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            _run_multi_sample_mode(args)
            output = fake_out.getvalue()
        self.assertIn("Erro", output)

    def test_multi_sample_mode_only_reference(self):
        """Test _run_multi_sample_mode with only one sequence shows warning."""
        from main import _run_multi_sample_mode, parse_args
        path = self.create_temp_fasta(">ref\nACTG")
        args = parse_args(["--input", path])
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            _run_multi_sample_mode(args)
            output = fake_out.getvalue()
        self.assertIn("Aviso", output)

    def test_multi_sample_mode_full_flow(self):
        """Test _run_multi_sample_mode with reference + samples."""
        from main import _run_multi_sample_mode, parse_args
        path = self.create_temp_fasta(">ref\nACTG\n>s1\nACTT\n>s2\nACTG")
        args = parse_args(["--input", path])
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            _run_multi_sample_mode(args)
            output = fake_out.getvalue()
        self.assertIn("s1", output)
        self.assertIn("s2", output)
        # Clean up generated report files
        for f in ["snps_report_s1.txt"]:
            if os.path.exists(f):
                os.remove(f)


if __name__ == "__main__":
    unittest.main()
