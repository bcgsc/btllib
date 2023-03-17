import os
import btllib
import unittest


class SeqReaderTests(unittest.TestCase):

    def setUp(self):
        self.base_dir = os.path.dirname(__file__)

    def test_indexlr_single_seq(self):
        minimizers = [9638000835411166970, 16169464608842636080, 15783021188953468643]
        seq_id = "1"
        path = os.path.join(self.base_dir, 'single_seq.fa')
        with btllib.Indexlr(path, 10, 5, btllib.IndexlrFlag.LONG_MODE) as indexlr:
            for record in indexlr:
                self.assertEqual(seq_id, record.id)
                self.assertEqual(minimizers, [mx.out_hash for mx in record.minimizers])
