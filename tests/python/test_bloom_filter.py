import btllib
import unittest

class BloomFilterTests(unittest.TestCase):
    
    def test_seq_insertion(self):
        seq, k = 'AGTCATCGACTGATGC', 5
        bf = btllib.KmerBloomFilter(1024, 3, k)
        bf.insert(seq)
        self.assertEqual(bf.contains('TCATC'), 1)
