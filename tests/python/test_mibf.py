import btllib
import unittest


class MIBloomFilterTests(unittest.TestCase):

    def test_mibf_insertion(self):
        mibf = btllib.MIBloomFilter8(1024 * 1024, 3)
        mibf.insert_bv([1, 10, 100])
        mibf.insert_bv([100, 200, 300])
        mibf.complete_bv_insertion()

        self.assertTrue(mibf.bv_contains([1, 10, 100]))
        self.assertTrue(mibf.bv_contains([100, 200, 300]))
        self.assertFalse(mibf.bv_contains([1, 20, 100]))

        inserted_id = 12
        mibf.insert_id([1, 10, 100], inserted_id)
        for query_id in mibf.get_id([1, 10, 100]):
            self.assertEqual(query_id, inserted_id)
