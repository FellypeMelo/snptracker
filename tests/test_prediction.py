"""Tests for prediction.py — Grantham Score functional prediction module."""

import unittest
from prediction import (
    translate_dna_to_protein,
    get_amino_acid_change,
    grantham_score,
    grantham_prediction,
    predict_functional_impact,
)


class TestTranslateDnaToProtein(unittest.TestCase):

    def test_frame_positive_1(self):
        """Frame +1: translates from position 1."""
        result = translate_dna_to_protein("ATGGTTGCT", frame=1)
        self.assertEqual(result, ["Met", "Val", "Ala"])

    def test_frame_positive_2(self):
        """Frame +2: first base ignored, codons start at position 2."""
        result = translate_dna_to_protein("AATGGTTGCT", frame=2)
        self.assertEqual(result, ["Met", "Val", "Ala"])

    def test_frame_negative_1(self):
        """Frame -1: sequence is reverse-complemented then translated."""
        # rev_comp("TTACCAT") = "ATGGTAA" -> ATG(Met)+GTA(Val)
        result = translate_dna_to_protein("TTACCAT", frame=-1)
        self.assertEqual(result, ["Met", "Val"])

    def test_stop_codon_terminates_translation(self):
        """STOP codon is represented as '*' and ends the list."""
        # ATG(Met) TAA(STOP) GCT(never reached)
        result = translate_dna_to_protein("ATGTAAGCT", frame=1)
        self.assertEqual(result, ["Met", "*"])

    def test_incomplete_codon_at_end_is_ignored(self):
        """Trailing bases that do not form a full codon are ignored."""
        # ATG(Met) GTT(Val) + one leftover base
        result = translate_dna_to_protein("ATGGTTA", frame=1)
        self.assertEqual(result, ["Met", "Val"])


class TestGetAminoAcidChange(unittest.TestCase):
    # Reference: "ATGGTTGCT" -> Met(1) Val(2) Ala(3)
    # NON_SYNONYMOUS: pos 4, G->A -> GTT(Val)->ATT(Ile)
    _REF = "ATGGTTGCT"
    _NON_SYN_SNP = {
        "position": 4,
        "reference": "G",
        "alternate": "A",
        "type": "TRANSITION",
        "annotation": "NON_SYNONYMOUS",
        "context": "G[G>A]T",
    }

    def test_non_synonymous_returns_tuple(self):
        """NON_SYNONYMOUS SNP returns (ref_aa, aa_pos, alt_aa)."""
        result = get_amino_acid_change(self._NON_SYN_SNP, self._REF, frame=1)
        self.assertEqual(result, ("Val", 2, "Ile"))

    def test_synonymous_returns_none(self):
        """SYNONYMOUS SNP returns None — no amino acid change to score."""
        snp = {
            "position": 9, "reference": "T", "alternate": "C",
            "type": "TRANSITION", "annotation": "SYNONYMOUS",
            "context": "C[T>C]_",
        }
        self.assertIsNone(get_amino_acid_change(snp, self._REF, frame=1))

    def test_nonsense_returns_none(self):
        """NONSENSE SNP returns None — introduces stop, no substitution score."""
        snp = {
            "position": 4, "reference": "G", "alternate": "T",
            "type": "TRANSVERSION", "annotation": "NONSENSE",
            "context": "G[G>T]T",
        }
        self.assertIsNone(get_amino_acid_change(snp, self._REF, frame=1))

    def test_non_coding_returns_none(self):
        """NON_CODING SNP returns None."""
        snp = {
            "position": 1, "reference": "A", "alternate": "C",
            "type": "TRANSVERSION", "annotation": "NON_CODING",
            "context": "_[A>C]T",
        }
        self.assertIsNone(get_amino_acid_change(snp, self._REF, frame=1))

    def test_indel_returns_none(self):
        """DELETION (INDEL) returns None — no substitution to score."""
        snp = {
            "position": 7,
            "reference": "-",
            "alternate": "A",
            "type": "INSERTION",
            "annotation": "NON_CODING",
        }
        self.assertIsNone(get_amino_acid_change(snp, self._REF, frame=1))


class TestGranthamScore(unittest.TestCase):

    def test_conservative_pair_val_ile(self):
        """Val->Ile is a well-known conservative substitution (score 29)."""
        self.assertEqual(grantham_score("Val", "Ile"), 29)

    def test_radical_pair_cys_trp(self):
        """Cys->Trp is a radical substitution (score > 150)."""
        score = grantham_score("Cys", "Trp")
        self.assertGreater(score, 150)

    def test_identical_amino_acids_score_zero(self):
        """Same amino acid on both sides yields score 0."""
        self.assertEqual(grantham_score("Val", "Val"), 0)

    def test_symmetry(self):
        """Grantham score is symmetric: d(A,B) == d(B,A)."""
        self.assertEqual(grantham_score("Ala", "Glu"), grantham_score("Glu", "Ala"))

    def test_stop_codon_returns_none(self):
        """STOP is not a standard amino acid — score returns None."""
        self.assertIsNone(grantham_score("Val", "STOP"))
        self.assertIsNone(grantham_score("STOP", "Val"))


class TestGranthamPrediction(unittest.TestCase):

    def test_score_zero_is_conservative(self):
        self.assertEqual(grantham_prediction(0), "CONSERVATIVE")

    def test_score_50_is_conservative(self):
        self.assertEqual(grantham_prediction(50), "CONSERVATIVE")

    def test_score_51_is_moderate(self):
        self.assertEqual(grantham_prediction(51), "MODERATE")

    def test_score_100_is_moderate(self):
        self.assertEqual(grantham_prediction(100), "MODERATE")

    def test_score_101_is_radical(self):
        self.assertEqual(grantham_prediction(101), "RADICAL")

    def test_very_high_score_is_radical(self):
        self.assertEqual(grantham_prediction(215), "RADICAL")


class TestPredictFunctionalImpact(unittest.TestCase):
    # Reference: "ATGGTTGCT" -> Met(1) Val(2) Ala(3)
    _REF = "ATGGTTGCT"

    def test_non_synonymous_returns_grantham_keys(self):
        """NON_SYNONYMOUS SNP returns dict with grantham_score and prediction."""
        # pos 4, G->A: GTT(Val)->ATT(Ile), score=29, CONSERVATIVE
        snp = {
            "position": 4, "reference": "G", "alternate": "A",
            "type": "TRANSITION", "annotation": "NON_SYNONYMOUS",
            "context": "G[G>A]T",
        }
        result = predict_functional_impact(snp, self._REF, frame=1)
        self.assertEqual(result["grantham_score"], 29)
        self.assertEqual(result["grantham_prediction"], "CONSERVATIVE")

    def test_synonymous_returns_empty_dict(self):
        """SYNONYMOUS SNP returns empty dict — no prediction to compute."""
        snp = {
            "position": 9, "reference": "T", "alternate": "C",
            "type": "TRANSITION", "annotation": "SYNONYMOUS",
            "context": "C[T>C]_",
        }
        self.assertEqual(predict_functional_impact(snp, self._REF, frame=1), {})

    def test_non_coding_returns_empty_dict(self):
        """NON_CODING SNP returns empty dict."""
        snp = {
            "position": 1, "reference": "A", "alternate": "C",
            "type": "TRANSVERSION", "annotation": "NON_CODING",
            "context": "_[A>C]T",
        }
        self.assertEqual(predict_functional_impact(snp, self._REF, frame=1), {})

    def test_radical_substitution_cys_trp(self):
        """Cys->Trp yields RADICAL prediction."""
        # "ATGTGTGCT" -> Met(1) Cys(2) Ala(3)
        # pos 6, T->G: TGT(Cys)->TGG(Trp)
        ref = "ATGTGTGCT"
        snp = {
            "position": 6, "reference": "T", "alternate": "G",
            "type": "TRANSVERSION", "annotation": "NON_SYNONYMOUS",
            "context": "T[T>G]G",
        }
        result = predict_functional_impact(snp, ref, frame=1)
        self.assertEqual(result["grantham_prediction"], "RADICAL")
        self.assertGreater(result["grantham_score"], 150)


if __name__ == "__main__":
    unittest.main()
