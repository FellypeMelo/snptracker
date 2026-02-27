import unittest
import os
from unittest.mock import patch
import io
from main import detect_snps, classify_mutation, generate_snp_file, print_snp_report

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

if __name__ == "__main__":
    unittest.main()
