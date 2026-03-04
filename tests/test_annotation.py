import unittest
from annotation import annotate_snp, get_codon, translate_codon


class TestTranslateCodon(unittest.TestCase):

    def test_translate_codon_alanine(self):
        """GCT encodes Alanine."""
        self.assertEqual(translate_codon("GCT"), "Ala")

    def test_translate_codon_alanine_synonymous_variant(self):
        """GCC also encodes Alanine (synonymous)."""
        self.assertEqual(translate_codon("GCC"), "Ala")

    def test_translate_codon_stop_taa(self):
        """TAA is a stop codon."""
        self.assertEqual(translate_codon("TAA"), "STOP")

    def test_translate_codon_stop_tag(self):
        """TAG is a stop codon."""
        self.assertEqual(translate_codon("TAG"), "STOP")

    def test_translate_codon_stop_tga(self):
        """TGA is a stop codon."""
        self.assertEqual(translate_codon("TGA"), "STOP")

    def test_translate_codon_methionine(self):
        """ATG encodes Methionine (start codon)."""
        self.assertEqual(translate_codon("ATG"), "Met")

    def test_translate_codon_invalid_too_short(self):
        """String shorter than 3 chars raises ValueError."""
        with self.assertRaises(ValueError):
            translate_codon("GC")

    def test_translate_codon_invalid_chars(self):
        """Codon with non-ACGT characters raises ValueError."""
        with self.assertRaises(ValueError):
            translate_codon("GCN")

    def test_translate_codon_lowercase_raises(self):
        """Lowercase codon raises ValueError (expects uppercase)."""
        with self.assertRaises(ValueError):
            translate_codon("gct")


class TestGetCodon(unittest.TestCase):

    def test_get_codon_first_triplet(self):
        """Position 1 in 'GCTAAA' returns first codon 'GCT'."""
        self.assertEqual(get_codon("GCTAAA", 1), "GCT")

    def test_get_codon_second_triplet_start(self):
        """Position 4 in 'GCTAAA' returns second codon 'AAA'."""
        self.assertEqual(get_codon("GCTAAA", 4), "AAA")

    def test_get_codon_second_triplet_middle(self):
        """Position 5 in 'GCTAAA' also returns second codon 'AAA'."""
        self.assertEqual(get_codon("GCTAAA", 5), "AAA")

    def test_get_codon_second_triplet_end(self):
        """Position 6 in 'GCTAAA' also returns second codon 'AAA'."""
        self.assertEqual(get_codon("GCTAAA", 6), "AAA")

    def test_get_codon_incomplete_trailing_codon(self):
        """Position 7 in 'GCTAAAG' falls in incomplete codon, returns ''."""
        self.assertEqual(get_codon("GCTAAAG", 7), "")

    def test_get_codon_sequence_shorter_than_3(self):
        """Position 1 in 'AC' (< 3 bases) returns ''."""
        self.assertEqual(get_codon("AC", 1), "")

    def test_get_codon_position_zero_raises(self):
        """Position 0 raises ValueError."""
        with self.assertRaises(ValueError):
            get_codon("GCTAAA", 0)

    def test_get_codon_position_exceeds_length_raises(self):
        """Position beyond sequence length raises ValueError."""
        with self.assertRaises(ValueError):
            get_codon("GCT", 10)


class TestAnnotateSnp(unittest.TestCase):

    def test_annotate_snp_synonymous(self):
        """GCT→GCC: both encode Ala → SYNONYMOUS."""
        ref = "GCTAAA"
        alt = "GCCAAA"
        self.assertEqual(annotate_snp(3, ref, alt), "SYNONYMOUS")

    def test_annotate_snp_non_synonymous(self):
        """GCT→GTT: Ala→Val → NON_SYNONYMOUS."""
        ref = "GCTAAA"
        alt = "GTTAAA"
        self.assertEqual(annotate_snp(2, ref, alt), "NON_SYNONYMOUS")

    def test_annotate_snp_nonsense(self):
        """GCT→TAA: Ala→STOP → NONSENSE."""
        ref = "GCTAAA"
        alt = "TAAAAA"
        self.assertEqual(annotate_snp(1, ref, alt), "NONSENSE")

    def test_annotate_snp_non_coding_incomplete_codon(self):
        """SNP in trailing partial codon returns NON_CODING."""
        ref = "GCTAAAG"
        alt = "GCTAAAT"
        self.assertEqual(annotate_snp(7, ref, alt), "NON_CODING")

    def test_annotate_snp_non_coding_sequence_too_short(self):
        """Sequence shorter than 3 bases returns NON_CODING."""
        self.assertEqual(annotate_snp(1, "AC", "AT"), "NON_CODING")

    def test_annotate_snp_all_stop_codons_are_nonsense(self):
        """Any substitution producing TAG or TGA is also NONSENSE."""
        # GCT→TAG
        ref_tag = "GCTAAA"
        alt_tag = "TAGAAA"
        self.assertEqual(annotate_snp(1, ref_tag, alt_tag), "NONSENSE")

        # GCT→TGA
        alt_tga = "TGAAAA"
        self.assertEqual(annotate_snp(1, ref_tag, alt_tga), "NONSENSE")

    def test_annotate_snp_second_codon_non_synonymous(self):
        """SNP in the second codon is annotated correctly."""
        # AAA→AAG: both Lys → SYNONYMOUS
        ref = "GCTAAA"
        alt = "GCTAAG"
        self.assertEqual(annotate_snp(6, ref, alt), "SYNONYMOUS")


if __name__ == "__main__":
    unittest.main()
