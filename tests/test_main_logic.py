import unittest
import os
from unittest.mock import patch
import io
from main import (
    detect_snps,
    classify_mutation,
    generate_snp_file,
    print_snp_report,
    run_multi_sample,
    get_trinucleotide_context,
    parse_cds_regions,
)

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

    def test_detect_snps_has_context_key(self):
        """SNP dicts must contain 'context' key."""
        snps = detect_snps("TACTG", "TACGG")
        snp = snps[0]
        self.assertIn('context', snp)

    def test_detect_snps_indel_has_no_context(self):
        """INSERTION and DELETION dicts must NOT contain 'context' key."""
        ins_snps = detect_snps("ACTG", "ACTGA")
        self.assertNotIn('context', ins_snps[0])

        del_snps = detect_snps("ACTG", "ACT")
        self.assertNotIn('context', del_snps[0])

    def test_detect_snps_context_value_correct(self):
        """Context value must match COSMIC format for a known SNP."""
        # ref="TACTG", smp="TACGG" → SNP at pos 4: T>G, flanks T and G... wait
        # ref="TACTG" pos 4 = T, smp pos 4 = G → prev=ref[2]='C', next=ref[4]='G'
        snps = detect_snps("TACTG", "TACGG")
        self.assertEqual(snps[0]['context'], "C[T>G]G")


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


class TestParseCdsRegions(unittest.TestCase):

    def test_single_region(self):
        """'1-90' parses to [(1, 90)]."""
        self.assertEqual(parse_cds_regions("1-90"), [(1, 90)])

    def test_multiple_regions(self):
        """'1-90,100-150' parses to [(1, 90), (100, 150)]."""
        self.assertEqual(
            parse_cds_regions("1-90,100-150"),
            [(1, 90), (100, 150)]
        )

    def test_invalid_format_raises(self):
        """Non-numeric input raises ValueError."""
        with self.assertRaises(ValueError):
            parse_cds_regions("abc")

    def test_start_greater_than_end_raises(self):
        """start > end raises ValueError."""
        with self.assertRaises(ValueError):
            parse_cds_regions("90-10")


class TestDetectSnpsWithCds(unittest.TestCase):

    def test_snp_outside_cds_returns_non_coding(self):
        """SNP outside CDS region gets annotation NON_CODING."""
        # ref="GCTGCT", pos 4 is G→T; region only covers (1, 3)
        snps = detect_snps("GCTGCT", "GCTTCT", cds_regions=[(1, 3)])
        self.assertEqual(snps[0]["annotation"], "NON_CODING")

    def test_snp_inside_cds_returns_functional_annotation(self):
        """SNP inside CDS region gets functional annotation."""
        # ref="GCTAAA", pos 2 → GCT→GTT → NON_SYNONYMOUS
        snps = detect_snps("GCTAAA", "GTTAAA", cds_regions=[(1, 6)])
        self.assertNotEqual(snps[0]["annotation"], "NON_CODING")

    def test_no_cds_regions_behaves_as_before(self):
        """detect_snps without cds_regions is fully backwards-compatible."""
        snps_old = detect_snps("GCTAAA", "GTTAAA")
        snps_new = detect_snps("GCTAAA", "GTTAAA", cds_regions=None)
        self.assertEqual(snps_old[0]["annotation"], snps_new[0]["annotation"])


class TestGetTrinucleotideContext(unittest.TestCase):

    def test_middle_position(self):
        """SNP in the middle returns correct flanking bases."""
        self.assertEqual(
            get_trinucleotide_context("TACTG", 3, "C", "G"),
            "A[C>G]T"
        )

    def test_left_edge(self):
        """SNP at position 1 uses '_' as left sentinel."""
        self.assertEqual(
            get_trinucleotide_context("ACTG", 1, "A", "G"),
            "_[A>G]C"
        )

    def test_right_edge(self):
        """SNP at last position uses '_' as right sentinel."""
        self.assertEqual(
            get_trinucleotide_context("ACTG", 4, "G", "T"),
            "T[G>T]_"
        )

    def test_single_base_sequence(self):
        """SNP in a single-base sequence uses '_' on both sides."""
        self.assertEqual(
            get_trinucleotide_context("A", 1, "A", "G"),
            "_[A>G]_"
        )

    def test_position_two(self):
        """SNP at position 2 returns correct flanking bases."""
        self.assertEqual(
            get_trinucleotide_context("ACTG", 2, "C", "T"),
            "A[C>T]T"
        )


class TestDetectSnpsFrame(unittest.TestCase):
    """Tests for detect_snps() with explicit reading frame."""

    def test_frame1_default_backwards_compatible(self):
        """detect_snps with no frame arg uses frame +1 — same as before."""
        ref = "GCTAAA"
        smp = "GCCAAA"  # pos3 T→C; frame+1 codon GCT→GCC (both Ala): SYNONYMOUS
        snps = detect_snps(ref, smp)
        self.assertEqual(len(snps), 1)
        self.assertEqual(snps[0]["annotation"], "SYNONYMOUS")

    def test_frame2_annotates_correctly(self):
        """detect_snps with frame=2 annotates using frame +2 codon."""
        # GCTAAACGT; frame+2: CTA|AAC|GT
        # SNP at pos5: A→C → second codon AAC→CAC (Asn→His): NON_SYNONYMOUS
        ref = "GCTAAACGT"
        smp = "GCTACACGT"  # pos5 A→C only
        snps = detect_snps(ref, smp, frame=2)
        self.assertEqual(len(snps), 1)
        self.assertEqual(snps[0]["position"], 5)
        self.assertEqual(snps[0]["annotation"], "NON_SYNONYMOUS")

    def test_frame2_position1_is_non_coding(self):
        """SNP at pos 1 is NON_CODING in frame +2 (before frame start)."""
        ref = "GCTAAACGT"
        smp = "TCTAAACGT"  # pos1 G→T
        snps = detect_snps(ref, smp, frame=2)
        self.assertEqual(snps[0]["annotation"], "NON_CODING")

    def test_frame3_positions_1_2_non_coding(self):
        """SNPs at positions 1 and 2 are NON_CODING in frame +3."""
        ref = "GCTAAACGT"
        smp = "TATCAACGT"  # pos1 G→T, pos2 C→A
        snps = detect_snps(ref, smp, frame=3)
        positions_annotated = {s["position"]: s["annotation"] for s in snps}
        self.assertEqual(positions_annotated[1], "NON_CODING")
        self.assertEqual(positions_annotated[2], "NON_CODING")

    def test_frame_minus1_annotation(self):
        """detect_snps with frame=-1 annotates via reverse complement."""
        ref = "GCTAAACGT"
        smp = "GCTAAACGA"  # pos9 T→A
        snps = detect_snps(ref, smp, frame=-1)
        self.assertEqual(len(snps), 1)
        self.assertEqual(snps[0]["annotation"], "NON_SYNONYMOUS")

    def test_frame_with_cds_regions(self):
        """frame param is combined correctly with cds_regions."""
        ref = "GCTAAACGT"
        smp = "GCTACACGT"  # pos5 A→C only
        # CDS covers the whole sequence, frame=2
        snps = detect_snps(ref, smp, cds_regions=[(1, 9)], frame=2)
        self.assertEqual(snps[0]["annotation"], "NON_SYNONYMOUS")

    def test_frame_invalid_raises(self):
        """detect_snps propagates ValueError for invalid frame."""
        with self.assertRaises(ValueError):
            detect_snps("GCTAAA", "GCCAAA", frame=0)


if __name__ == "__main__":
    unittest.main()
