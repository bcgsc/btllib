import os
import unittest
import btllib
import random
import string

class TestCountingBloomFilter8(unittest.TestCase):

    def setUp(self):
        self.base_dir = os.path.dirname(__file__)

    def test_insert_and_contains(self):
        cbf = btllib.CountingBloomFilter8(1024 * 1024, 3)
        element = ([1, 10, 100])
        cbf.insert(element)
        cbf.insert(element)
        self.assertEqual(cbf.contains(element), 2)

    def test_absent_element(self):
        cbf = btllib.CountingBloomFilter8(1024 * 1024, 3)
        element = ([1, 10, 100])
        cbf.insert(element)
        cbf.insert(element)

        absent_element = ([1, 20, 100])
        self.assertEqual(cbf.contains(absent_element), 0)

    def test_save_load(self):
        cbf = btllib.CountingBloomFilter8(1024 * 1024, 3)
        element = ([1, 10, 100])
        cbf.insert(element)
        cbf.insert(element)
        
        file_path = os.path.join(self.base_dir, "test.cbf")
        cbf.save(file_path)
        try:
            loaded_cbf = btllib.CountingBloomFilter8(file_path)
            element = ([1, 10, 100])
            self.assertEqual(loaded_cbf.contains(element), 2)
        finally:
            os.remove(file_path)
        
