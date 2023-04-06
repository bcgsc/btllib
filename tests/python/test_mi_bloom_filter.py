import os
import random
import string
import btllib
import unittest


class MIBloomFilterTests(unittest.TestCase):

    def setUp(self):
        self.base_dir = os.path.dirname(__file__)

    def test_mibloomfilter_bv_and_id_insertion(self):
        mi_bf_1 = btllib.MIBloomFilter8(1024 * 1024, 3, "ntHash")
        mi_bf_1.insert_bv([1, 10, 100])
        mi_bf_1.insert_bv([100, 200, 300])
        mi_bf_1.complete_bv_insertion()

        self.assertTrue(mi_bf_1.bv_contains([1, 10, 100]))
        self.assertTrue(mi_bf_1.bv_contains([100, 200, 300]))
        self.assertFalse(mi_bf_1.bv_contains([1, 20, 100]))

        ID_1 = 12
        mi_bf_1.insert_id([1, 10, 100], ID_1)

        results_1 = mi_bf_1.get_id([1, 10, 100])
        for id in results_1:
            self.assertEqual(id, ID_1)

        print("multi-indexed BloomFilter ID count test")
        print("Testing ID counting")
        include_saturated = True

        ID_count = mi_bf_1.get_id_occurence_count(include_saturated)
        self.assertEqual(ID_count[ID_1], 3)

    def test_mibloomfilter_random_sampling(self):
        print("Testing multi-indexed BloomFilter random sampling")
        
        random_dna = "GGTAGACACACGTCCACCCCGCTGCTCTGTGACAGGGACTAAAGAGGCGAAGATTATCGTGTGTGCCCCGTTATGGTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCGAATTCCGAGCGATTAAGCGTGACAGTCCCAGCGAACCCACAAAACGTGATCGCAGTCCATGCGATCATACGCAAGAAGGAAGGTCCCCATACACCGACGCACCAGTTTACACGCCGTATGCATAAACGAGCTGCACAAACGAGAGTGCTTGAACTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAGTTGCCTTGCCTAGATGCAATGTCGGACGTATTACTTTTGCCTCAACGGCTCCTGCTTTCGCTGAAACCCAAGACAGGCAACAGTAACCGCCTTTTGAAGGCGAGTCCTTCGTCTGTGACTAACTGTGCCAAATCGTCTTCCAAACTCCTAATCCAGTTTAACTCACCAAATTATAGCCATACAGACCCTAATTTCATATCATATCACGCCATTAGCCTCTGCTAAAATTCTGTGCTCAAGGGTTTTGGTTCGCCCGAGTGATGTTGCCAATTAGGACCATCAAATGCACATGTTACAGGACTTCTTATAAATACTTTTTTCCTGGGGAGTAGCGGATCTTAATGGATGTTGCCAGCTGGTATGGAAGCTAATAGCGCCGGTGGGAGCGTAATCTGCCGTCTCCACCAACACAACGCTATCGGGTCATATTATAAGATTCCGCAATGGGGTTACTTATAGGTAGCCTTAACGATATCCGGAACTTGCGATGTACGTGCTATGCTTTAATACATACCTGGCCCAGTAGTTTTCCAATATGGGAACATCAATTGTACATCGGGCCGGGATAATCATGTCATCACGGAAGTAGCCGTAAGACAAATAATTCAAAAGAGATGTCGTTTTGCTAGTTCACGTGAAGGTGTCTCGCGCCACCTCTAAGTAAGTGGGCCGTCGAGA"
        dna_length = len(random_dna)
        expected_id_count = dna_length // 4
        tolerance = 0.1

        mi_bf = btllib.MIBloomFilter8(1024 * 1024, 1, "ntHash")
        nthash = btllib.NtHash(random_dna, 1, 15)
        while nthash.roll():
            mi_bf.insert_bv(nthash.hashes())

        mi_bf.complete_bv_insertion()

        ID_array = {0, 1, 2, 3}
        for id in ID_array:
            nthash = btllib.NtHash(random_dna, 1, 15)
            while nthash.roll():
                mi_bf.insert_id(nthash.hashes(), id)

        results_2 = [0] * 1
        total_counter = [0] * 4

        nthash = btllib.NtHash(random_dna, 1, 15)
        while nthash.roll():
            results_2 = mi_bf.get_id(nthash.hashes())
            for res in results_2:
                total_counter[res] += 1

        for count in total_counter:
            assert count < expected_id_count + (expected_id_count * tolerance)
            assert count > expected_id_count - (expected_id_count * tolerance)

        assert mi_bf.get_pop_saturated_cnt() == 0
