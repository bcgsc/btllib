import os
import random
import string
import btllib
import unittest

class MIBloomFilterTests(unittest.TestCase):

    def setUp(self):
        self.base_dir = os.path.dirname(__file__)
        
        ### mi_bf_1
        self.test_hashes_1 = [[1, 10, 100], [100, 200, 300]]
        ###########
        
        ### mi_bf_2
        self.random_dna_2 = "GGTAGACACACGTCCACCCCGCTGCTCTGTGACAGGGACTAAAGAGGCGAAGATTATCGTGTGTGCCCCGTTATGGTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCGAATTCCGAGCGATTAAGCGTGACAGTCCCAGCGAACCCACAAAACGTGATCGCAGTCCATGCGATCATACGCAAGAAGGAAGGTCCCCATACACCGACGCACCAGTTTACACGCCGTATGCATAAACGAGCTGCACAAACGAGAGTGCTTGAACTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAGTTGCCTTGCCTAGATGCAATGTCGGACGTATTACTTTTGCCTCAACGGCTCCTGCTTTCGCTGAAACCCAAGACAGGCAACAGTAACCGCCTTTTGAAGGCGAGTCCTTCGTCTGTGACTAACTGTGCCAAATCGTCTTCCAAACTCCTAATCCAGTTTAACTCACCAAATTATAGCCATACAGACCCTAATTTCATATCATATCACGCCATTAGCCTCTGCTAAAATTCTGTGCTCAAGGGTTTTGGTTCGCCCGAGTGATGTTGCCAATTAGGACCATCAAATGCACATGTTACAGGACTTCTTATAAATACTTTTTTCCTGGGGAGTAGCGGATCTTAATGGATGTTGCCAGCTGGTATGGAAGCTAATAGCGCCGGTGGGAGCGTAATCTGCCGTCTCCACCAACACAACGCTATCGGGTCATATTATAAGATTCCGCAATGGGGTTACTTATAGGTAGCCTTAACGATATCCGGAACTTGCGATGTACGTGCTATGCTTTAATACATACCTGGCCCAGTAGTTTTCCAATATGGGAACATCAATTGTACATCGGGCCGGGATAATCATGTCATCACGGAAGTAGCCGTAAGACAAATAATTCAAAAGAGATGTCGTTTTGCTAGTTCACGTGAAGGTGTCTCGCGCCACCTCTAAGTAAGTGGGCCGTCGAGA"
        self.id_array_2 = {0, 1, 2, 3}
        ###########
   
    def set_up_mi_bf_1(self):
        self.mi_bf_1 = btllib.MIBloomFilter8(1024 * 1024, 3, "ntHash")
        for h in self.test_hashes_1:
            self.mi_bf_1.insert_bv(h)
        self.mi_bf_1.complete_bv_insertion()

    def set_up_mi_bf_2(self):
        self.mi_bf_2 = btllib.MIBloomFilter8(1024 * 1024, 3, "ntHash")
        
        nthash = btllib.NtHash(self.random_dna_2, 3, 15)
        while nthash.roll():
            self.mi_bf_2.insert_bv(nthash.hashes())
        self.mi_bf_2.complete_bv_insertion()
        
        for id in self.id_array_2:
            nthash = btllib.NtHash(self.random_dna_2, 3, 15)
            while nthash.roll():
                self.mi_bf_2.insert_id(nthash.hashes(), id)
         
    def test_mibloomfilter_bv_insertion(self):
        
        self.set_up_mi_bf_1() 

        self.assertTrue(self.mi_bf_1.bv_contains([1, 10, 100]))
        self.assertTrue(self.mi_bf_1.bv_contains([100, 200, 300]))
        self.assertFalse(self.mi_bf_1.bv_contains([1, 20, 100]))

    def test_mibloomfilter_id_insertion(self):

        self.set_up_mi_bf_1()
        
        expected_id = 12
        self.mi_bf_1.insert_id(self.test_hashes_1[0], expected_id)

        for result_id in self.mi_bf_1.get_id(self.test_hashes_1[0]):
            self.assertEqual(result_id, expected_id)
            
    def test_mibloomfilter_id_occurence(self):
        self.set_up_mi_bf_1()
        
        expected_id = 12
        self.mi_bf_1.insert_id(self.test_hashes_1[0], expected_id)

        include_saturated = True
        self.assertEqual(len(self.test_hashes_1[0]), 
                         self.mi_bf_1.get_id_occurence_count(include_saturated)[expected_id])

    def test_mibloomfilter_random_sampling(self):
        
        self.set_up_mi_bf_2()

        id_occurences = self.mi_bf_2.get_id_occurence_count(False)[0:len(self.id_array_2)] ## function returns fixed size large array with zeroes
        
        total_id_occurences = sum(id_occurences)
        expected_id_occurence = total_id_occurences // len(self.id_array_2) 

        for occurence in id_occurences:
            self.assertAlmostEqual(occurence, 
                                   expected_id_occurence,
                                   delta = expected_id_occurence // 10)
            
    def test_mibloomfilter_saving(self):
        self.set_up_mi_bf_2()
        
        file_path = os.path.join(self.base_dir, "test.mibf")
        self.mi_bf_2.save(file_path)
       
        try:
            loaded_mi_bf = btllib.MIBloomFilter8(file_path)
        
            expected_id_occurences = self.mi_bf_2.get_id_occurence_count(True)
            id_occurences_loaded = loaded_mi_bf.get_id_occurence_count(True)
            
            for i in range(len(self.id_array_2)):
                self.assertEqual(id_occurences_loaded[i], expected_id_occurences[i])
        finally:
            os.remove(file_path)
            os.remove(file_path+".sdsl")
