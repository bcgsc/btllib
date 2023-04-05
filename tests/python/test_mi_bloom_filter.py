import os
import random
import string
import btllib
import unittest


class MIBloomFilterTests(unittest.TestCase):

    def setUp(self):
        self.base_dir = os.path.dirname(__file__)

    @staticmethod
    def swigPyObject_iter(swig_obj):
        i = 0
        while True:
            try:
                yield swig_obj[i]
                i += 1
            except IndexError:
                break

    def test_mibloomfilter_operations(self):
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
        results_1_list = self.swigPyObject_iter(results_1)
        for id in results_1_list:
            self.assertEqual(id, ID_1)
