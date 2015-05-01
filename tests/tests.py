import unittest

from AMAS import DNAAlignment


class TestAMAS(unittest.TestCase):
    def test_statistics(self):
        in_file = 'tests/fasta1.fas'
        in_format = 'fasta'
        data_type = 'dna'

        expected = ['10', '100', '1000', '1', '0.1', '2', '0.02', '1', '0.01']

        aln = DNAAlignment(in_file, in_format, data_type).summarize_alignment()
        result = aln[1:]
        self.assertEqual(expected, result)
