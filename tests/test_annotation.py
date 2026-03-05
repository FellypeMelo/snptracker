import unittest
from annotation import (
    annotate_snp,
    annotate_snp_with_regions,
    get_codon,
    is_in_cds,
    reverse_complement,
    translate_codon,
)


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


class TestIsInCds(unittest.TestCase):

    def test_position_inside_single_region(self):
        """Position within the only CDS region returns True."""
        self.assertTrue(is_in_cds(5, [(1, 10)]))

    def test_position_outside_single_region(self):
        """Position beyond the CDS region returns False."""
        self.assertFalse(is_in_cds(15, [(1, 10)]))

    def test_empty_cds_list(self):
        """No CDS regions defined — every position is NON_CODING."""
        self.assertFalse(is_in_cds(5, []))

    def test_position_on_exact_boundaries(self):
        """Start and end positions are inclusive."""
        self.assertTrue(is_in_cds(1, [(1, 10)]))
        self.assertTrue(is_in_cds(10, [(1, 10)]))

    def test_position_in_second_of_two_regions(self):
        """Position inside the second region returns True."""
        self.assertTrue(is_in_cds(80, [(10, 50), (70, 120)]))

    def test_position_between_two_regions(self):
        """Position in the gap between two regions returns False."""
        self.assertFalse(is_in_cds(60, [(10, 50), (70, 120)]))


class TestAnnotateSnpWithRegions(unittest.TestCase):

    def test_none_regions_behaves_like_annotate_snp(self):
        """cds_regions=None falls back to annotate_snp() behavior."""
        ref = "GCTAAA"
        alt = "GCCAAA"
        self.assertEqual(
            annotate_snp_with_regions(3, ref, alt, cds_regions=None),
            "SYNONYMOUS"
        )

    def test_position_outside_cds_is_non_coding(self):
        """SNP outside every CDS region returns NON_CODING."""
        ref = "GCTAAAGTT"
        alt = "GCTAATGTT"
        # pos 6 is outside region (1, 3)
        self.assertEqual(
            annotate_snp_with_regions(6, ref, alt, cds_regions=[(1, 3)]),
            "NON_CODING"
        )

    def test_position_inside_cds_annotated_normally(self):
        """SNP inside a CDS region gets full functional annotation."""
        ref = "GCTAAA"
        alt = "GTTAAA"
        # pos 2 → GCT→GTT → Ala→Val → NON_SYNONYMOUS
        self.assertEqual(
            annotate_snp_with_regions(2, ref, alt, cds_regions=[(1, 6)]),
            "NON_SYNONYMOUS"
        )

    def test_empty_cds_list_returns_non_coding(self):
        """cds_regions=[] forces every SNP to NON_CODING."""
        ref = "GCTAAA"
        alt = "GTTAAA"
        self.assertEqual(
            annotate_snp_with_regions(2, ref, alt, cds_regions=[]),
            "NON_CODING"
        )


class TestReverseComplement(unittest.TestCase):

    def test_single_base_a(self):
        """A complements to T and reverses to T."""
        self.assertEqual(reverse_complement("A"), "T")

    def test_single_base_c(self):
        self.assertEqual(reverse_complement("C"), "G")

    def test_single_base_g(self):
        self.assertEqual(reverse_complement("G"), "C")

    def test_single_base_t(self):
        self.assertEqual(reverse_complement("T"), "A")

    def test_complement_and_reverse(self):
        """ATCG → complement TAGC → reverse CGAT."""
        self.assertEqual(reverse_complement("ATCG"), "CGAT")

    def test_palindrome(self):
        """ATAT → complement TATA → reverse ATAT."""
        self.assertEqual(reverse_complement("ATAT"), "ATAT")

    def test_known_codon_atg(self):
        """ATG (Met) reverse complement is CAT (His)."""
        self.assertEqual(reverse_complement("ATG"), "CAT")

    def test_empty_string(self):
        self.assertEqual(reverse_complement(""), "")

    def test_invalid_base_raises(self):
        """Non-ACGT character in sequence raises ValueError."""
        with self.assertRaises(ValueError):
            reverse_complement("ACGN")


class TestGetCodonWithFrame(unittest.TestCase):
    # sequence: GCTAAACGT  (frame +1 codons: GCT|AAA|CGT)
    #                       (frame +2 codons: _CT|AAA|CGT  → first in-frame: pos2)
    #                       (frame +3 codons: __T|AAA|CGT  → first in-frame: pos3)
    SEQ = "GCTAAACGT"

    def test_frame2_position_before_start_is_empty(self):
        """Position 1 in frame +2 is before the frame start → ''."""
        self.assertEqual(get_codon(self.SEQ, 1, frame=2), "")

    def test_frame2_first_codon(self):
        """Position 2 in frame +2 starts first codon 'CTA'."""
        self.assertEqual(get_codon(self.SEQ, 2, frame=2), "CTA")

    def test_frame2_first_codon_last_base(self):
        """Position 4 in frame +2 still belongs to first codon 'CTA'."""
        self.assertEqual(get_codon(self.SEQ, 4, frame=2), "CTA")

    def test_frame2_second_codon(self):
        """Position 5 in frame +2 starts second codon 'AAC'."""
        self.assertEqual(get_codon(self.SEQ, 5, frame=2), "AAC")

    def test_frame3_positions_before_start_are_empty(self):
        """Positions 1 and 2 in frame +3 are before start → ''."""
        self.assertEqual(get_codon(self.SEQ, 1, frame=3), "")
        self.assertEqual(get_codon(self.SEQ, 2, frame=3), "")

    def test_frame3_first_codon(self):
        """Position 3 in frame +3 starts first codon 'TAA'."""
        self.assertEqual(get_codon(self.SEQ, 3, frame=3), "TAA")

    def test_frame3_second_codon(self):
        """Position 6 in frame +3 starts second codon 'ACG'."""
        self.assertEqual(get_codon(self.SEQ, 6, frame=3), "ACG")

    def test_frame1_unchanged(self):
        """Frame +1 explicit behaves like the default."""
        self.assertEqual(get_codon(self.SEQ, 1, frame=1), "GCT")
        self.assertEqual(get_codon(self.SEQ, 4, frame=1), "AAA")

    def test_invalid_frame_raises(self):
        """Frame value outside {1,2,3,-1,-2,-3} raises ValueError."""
        with self.assertRaises(ValueError):
            get_codon(self.SEQ, 1, frame=0)
        with self.assertRaises(ValueError):
            get_codon(self.SEQ, 1, frame=4)

    def test_frame_minus1_forward_position(self):
        """Frame -1: codon for position derived from reverse complement."""
        # SEQ = GCTAAACGT; rev_comp = ACGTTTAGC
        # frame -1, pos 9 → rev_pos=1 → first codon of ACGTTTAGC = ACG
        self.assertEqual(get_codon(self.SEQ, 9, frame=-1), "ACG")

    def test_frame_minus1_before_start_is_empty(self):
        """Frame -1: position whose rev_pos falls before frame start → ''."""
        # SEQ = GCTAAACGT (len=9); rev_comp = ACGTTTAGC
        # pos 9 → rev_pos=1 → ok
        # pos 1 → rev_pos=9 → in-frame (9-1)//3 = codon 2
        # pos 8 → rev_pos=2 → codon 0 start=0; ok
        # We need a position that triggers adjusted < 0 for frame -1.
        # frame -1 offset=0, so adjusted = rev_pos - 1 - 0 = rev_pos - 1
        # rev_pos >= 1 always, so adjusted >= 0 always for frame -1.
        # For frame -2, offset=1: pos whose rev_pos=1 → adjusted = -1 → ""
        seq = "GCTAAACGT"  # len=9
        # frame -2, rev_pos of pos=9 is 1 → adjusted = 1-1-1 = -1 → ""
        self.assertEqual(get_codon(seq, 9, frame=-2), "")

    def test_frame_minus2_second_position(self):
        """Frame -2: rev_pos=2 → first codon of rev_comp starting at offset 1."""
        # SEQ = GCTAAACGT; rev_comp = ACGTTTAGC
        # pos=8 → rev_pos=2; offset=1; adjusted=0; codon_start=1 → "CGT"
        self.assertEqual(get_codon(self.SEQ, 8, frame=-2), "CGT")


class TestAnnotateSnpFrame(unittest.TestCase):

    def test_frame2_synonymous(self):
        """Frame +2: SNP that keeps amino acid → SYNONYMOUS."""
        # frame +2, first codon starts at pos 2
        # ref seq: _CTAAA... → CTA=Leu; alt: _TTAAA → TTA=Leu → SYNONYMOUS
        ref = "GCTAAACGT"
        alt = "GTTAAACGT"  # pos 2: C→T; frame +2 codon CTA→TTA (both Leu)
        self.assertEqual(annotate_snp(2, ref, alt, frame=2), "SYNONYMOUS")

    def test_frame2_non_synonymous(self):
        """Frame +2: SNP that changes amino acid → NON_SYNONYMOUS."""
        # ref codon at pos2 in frame+2: CTA=Leu
        # change pos4 (last base of first codon): ref A→C → CTAvsCtC
        ref = "GCTAAACGT"
        alt = "GCTCAACGT"  # pos4 A→C; frame+2: CTA→CTC both Leu — let's pick differently
        # pos5: frame+2 second codon AAC, pos5=A; change to C → ACC
        # AAC=Asn, ACC=Thr → NON_SYNONYMOUS
        ref2 = "GCTAAACGT"
        alt2 = "GCTACCCGT"  # pos5 A→C; frame+2 codon AAC→ACC: Asn→Thr
        self.assertEqual(annotate_snp(5, ref2, alt2, frame=2), "NON_SYNONYMOUS")

    def test_frame2_nonsense(self):
        """Frame +2: SNP that introduces stop codon → NONSENSE."""
        # frame+2 first codon starts at pos2: CTA
        # Change pos2 C→T and pos3 T→A to get TAA? No, only single SNP.
        # Change pos4 (still in codon CTA): A→G → CTG=Leu, not stop.
        # Let's use a ref where the SNP creates a stop:
        # ref: _GCAAA → GCA=Ala; alt: _TAAAA → TAA=STOP → NONSENSE
        ref = "AGCAAACGT"   # frame+2: pos2 G starts → GCA|ACG|T
        alt = "ATAAAACGT"   # pos2 G→T → TCA=Ser, not stop
        # Better: craft codon that becomes stop
        # frame+2 codon1 = seq[1:4]; to get TAA need seq[1]T seq[2]A seq[3]A
        # ref: ATAAAA...; frame+2 codon1 = TAA = STOP already
        # Let's get NON_SYNONYMOUS → NONSENSE by single base:
        # ref codon1(frame+2) = TGG(Trp); alt changes last G→A → TGA=STOP
        ref3 = "ATGGCGTAA"  # frame+2: TGG=Trp, CGT=Arg
        alt3 = "ATGACGTAA"  # pos4 G→A; frame+2 codon1=TGA=STOP → NONSENSE
        self.assertEqual(annotate_snp(4, ref3, alt3, frame=2), "NONSENSE")

    def test_frame2_non_coding_before_frame_start(self):
        """Frame +2: SNP at position 1 is before frame start → NON_CODING."""
        ref = "GCTAAACGT"
        alt = "TCTAAACGT"  # pos1 G→T
        self.assertEqual(annotate_snp(1, ref, alt, frame=2), "NON_CODING")

    def test_frame3_non_coding_before_frame_start(self):
        """Frame +3: SNPs at positions 1 and 2 → NON_CODING."""
        ref = "GCTAAACGT"
        alt_pos1 = "TCTAAACGT"
        alt_pos2 = "GATAAACGT"
        self.assertEqual(annotate_snp(1, ref, alt_pos1, frame=3), "NON_CODING")
        self.assertEqual(annotate_snp(2, ref, alt_pos2, frame=3), "NON_CODING")

    def test_frame1_default_unchanged(self):
        """Explicit frame=1 gives same result as no frame argument."""
        ref = "GCTAAA"
        alt = "GCCAAA"
        self.assertEqual(annotate_snp(3, ref, alt), "SYNONYMOUS")
        self.assertEqual(annotate_snp(3, ref, alt, frame=1), "SYNONYMOUS")

    def test_frame_minus1_basic(self):
        """Frame -1: annotate SNP using reverse complement strand."""
        # SEQ = GCTAAACGT (len=9); rev_comp = ACGTTTAGC
        # frame -1, pos=9 → rev_pos=1 → codon ACG=Thr
        # Alternate: change pos=9 from T→A → alt=GCTAAACGA
        # rev_comp(alt) = TCGTTTAGC; rev_pos=1 → codon TCG=Ser
        # Thr → Ser: NON_SYNONYMOUS
        ref = "GCTAAACGT"
        alt = "GCTAAACGA"  # pos9 T→A
        self.assertEqual(annotate_snp(9, ref, alt, frame=-1), "NON_SYNONYMOUS")


class TestAnnotateSnpWithRegionsFrame(unittest.TestCase):

    def test_frame_passed_through_for_in_cds_position(self):
        """frame param is forwarded to annotate_snp when pos is in CDS."""
        ref = "GCTAAACGT"
        alt = "TCTAAACGT"  # pos1 G→T
        # frame=2 → pos1 before frame start → NON_CODING despite being in CDS
        result = annotate_snp_with_regions(
            1, ref, alt, cds_regions=[(1, 9)], frame=2
        )
        self.assertEqual(result, "NON_CODING")

    def test_frame_default_1_with_regions(self):
        """Default frame=1 with regions works as before."""
        ref = "GCTAAA"
        alt = "GTTAAA"
        self.assertEqual(
            annotate_snp_with_regions(2, ref, alt, cds_regions=[(1, 6)]),
            "NON_SYNONYMOUS"
        )

    def test_non_coding_position_ignores_frame(self):
        """Position outside CDS returns NON_CODING regardless of frame."""
        ref = "GCTAAAGTT"
        alt = "GCTAATGTT"
        self.assertEqual(
            annotate_snp_with_regions(
                6, ref, alt, cds_regions=[(1, 3)], frame=2
            ),
            "NON_CODING"
        )


if __name__ == "__main__":
    unittest.main()
