import unittest
import os
from fasta_parser import read_fasta

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

if __name__ == "__main__":
    unittest.main()
