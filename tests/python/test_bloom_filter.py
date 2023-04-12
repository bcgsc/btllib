import os
import random
import string
import btllib
import unittest

class BloomFilterTests(unittest.TestCase):
    
    def setUp(self):
        self.base_dir = os.path.dirname(__file__)
    
    def test_bloom_filter(self):
        bf = btllib.BloomFilter(1024 * 1024, 3, "ntHash")
        bf.insert([1, 10, 100])
        bf.insert([100, 200, 300])

        self.assertTrue(bf.contains([1, 10, 100]))
        self.assertTrue(bf.contains([100, 200, 300]))
        self.assertFalse(bf.contains([1, 20, 100]))

    def test_kmer_bloom_filter(self):
        seq = "CACTATCGACGATCATTCGAGCATCAGCGACTG"
        seq2 = "GTAGTACGATCAGCGACTATCGAGCTACGAGCA"

        kmer_bf = btllib.KmerBloomFilter(1024 * 1024, 4, len(seq) // 2)
        kmer_bf.insert(seq)
        self.assertEqual(kmer_bf.contains(seq), len(seq) - len(seq) // 2 + 1)
        self.assertLessEqual(kmer_bf.contains(seq2), 1)

    def test_seed_bloom_filter(self):
        seed1 = "000001111111111111111111111111111"
        seed2 = "111111111111111111111111111100000"
        seq = "CACTATCGACGATCATTCGAGCATCAGCGACTG"
        snp_seq1 = "AACTATCGACGATCATTCGAGCATCAGCGACTG"
        snp_seq2 = "CACTATCGACGATCATTCGAGCATCAGCGACTA"

        seed_bf = btllib.SeedBloomFilter(1024 * 1024, len(seq), [seed1, seed2], 4)
        seed_bf.insert(seq)
        hit_seeds = seed_bf.contains(seq)

        self.assertIn(0, hit_seeds[0])
        self.assertIn(1, hit_seeds[0])

        hit_seeds = seed_bf.contains(snp_seq1)
        self.assertIn(0, hit_seeds[0])
        self.assertNotIn(1, hit_seeds[0])

        hit_seeds = seed_bf.contains(snp_seq2)
        self.assertNotIn(0, hit_seeds[0])
        self.assertIn(1, hit_seeds[0])
        
    def test_bloom_filter_save_load(self):
        bf = btllib.BloomFilter(1024 * 1024, 3, "ntHash")
        bf.insert([1, 10, 100])
        bf.insert([100, 200, 300])

        file_path = os.path.join(self.base_dir, "test.bf")
        bf.save(file_path)

        try:
            bf2 = btllib.BloomFilter(file_path)

            self.assertEqual(bf2.get_hash_fn(), "ntHash")

            self.assertTrue(bf2.contains([1, 10, 100]))
            self.assertTrue(bf2.contains([100, 200, 300]))
            self.assertFalse(bf2.contains([1, 20, 100]))
        finally:
            os.remove(file_path)

