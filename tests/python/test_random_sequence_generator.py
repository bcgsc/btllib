import os
import btllib
import re
import unittest


class SeqReaderTests(unittest.TestCase):

    def setUp(self):
        self.base_dir = os.path.dirname(__file__)

    def test_random_sequence_dna(self):
        dna_re = re.compile(r"^[ACTG]*$")
        seq_generator = btllib.RandomSequenceGenerator(btllib.RandomSequenceGenerator.SequenceType_DNA)
        random_seq = seq_generator.generate(100)
        self.assertEqual(len(random_seq), 100)
        self.assertTrue(re.search(dna_re, random_seq))

    def test_random_sequence_rna(self):
        rna_re = re.compile(r"^[ACUG]*$")
        seq_generator = btllib.RandomSequenceGenerator(btllib.RandomSequenceGenerator.SequenceType_RNA)
        random_seq = seq_generator.generate(200)
        self.assertEqual(len(random_seq), 200)
        self.assertTrue(re.search(rna_re, random_seq))
