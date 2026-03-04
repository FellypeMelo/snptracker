import unittest
import os
from fasta_parser import read_fasta, read_all_sequences

class TestFastaParser(unittest.TestCase):
    def setUp(self):
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            if os.path.exists(f):
                os.remove(f)

    def create_temp_fasta(self, content):
        file_path = f"temp_{len(self.temp_files)}.fasta"
        with open(file_path, 'w') as f:
            f.write(content)
        self.temp_files.append(file_path)
        return file_path

    def test_read_fasta_single_line(self):
        """Test reading a simple single-line FASTA file."""
        path = self.create_temp_fasta(">seq1\nACTG")
        self.assertEqual(read_fasta(path), "ACTG")

    def test_read_fasta_multi_line(self):
        """Test reading a FASTA file with a sequence spread over multiple lines."""
        path = self.create_temp_fasta(">seq1\nACTG\nTAGC\n")
        self.assertEqual(read_fasta(path), "ACTGTAGC")

    def test_read_fasta_multiple_sequences(self):
        """Test that only the first sequence is returned from a multi-sequence file."""
        path = self.create_temp_fasta(">seq1\nACTG\n>seq2\nTAGC")
        self.assertEqual(read_fasta(path), "ACTG")

    def test_read_fasta_file_not_found(self):
        """Test that FileNotFoundError is raised for missing files."""
        with self.assertRaises(FileNotFoundError):
            read_fasta("non_existent.fasta")

    # --- Testes para read_all_sequences ---

    def test_read_all_sequences_multi(self):
        """Test reading multiple sequences returns all of them."""
        path = self.create_temp_fasta(">reference\nACTG\n>sample1\nACTT\n>sample2\nGCTG")
        result = read_all_sequences(path)
        self.assertEqual(len(result), 3)
        self.assertEqual(result[0], ("reference", "ACTG"))
        self.assertEqual(result[1], ("sample1", "ACTT"))
        self.assertEqual(result[2], ("sample2", "GCTG"))

    def test_read_all_sequences_single(self):
        """Test that a single-sequence file returns a list with one item."""
        path = self.create_temp_fasta(">ref\nACTG")
        result = read_all_sequences(path)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], ("ref", "ACTG"))

    def test_read_all_sequences_multiline(self):
        """Test that multi-line sequences are concatenated correctly."""
        path = self.create_temp_fasta(">ref\nACTG\nTAGC\n>smp\nGGGG")
        result = read_all_sequences(path)
        self.assertEqual(result[0], ("ref", "ACTGTAGC"))
        self.assertEqual(result[1], ("smp", "GGGG"))

    def test_read_all_sequences_empty_file(self):
        """Test that an empty file returns an empty list."""
        path = self.create_temp_fasta("")
        result = read_all_sequences(path)
        self.assertEqual(result, [])

    def test_read_all_sequences_file_not_found(self):
        """Test that FileNotFoundError is raised for missing files."""
        with self.assertRaises(FileNotFoundError):
            read_all_sequences("non_existent.fasta")

    def test_read_all_sequences_header_full_text(self):
        """Test that the full header text (after >) is captured."""
        path = self.create_temp_fasta(">reference Sample 1 of 3\nACTG")
        result = read_all_sequences(path)
        self.assertEqual(result[0][0], "reference Sample 1 of 3")

if __name__ == "__main__":
    unittest.main()
