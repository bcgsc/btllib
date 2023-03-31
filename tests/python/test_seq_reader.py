import os
import btllib
import unittest


class SeqReaderTests(unittest.TestCase):

    def setUp(self):
        self.base_dir = os.path.dirname(__file__)

    def test_seq_reader_single_seq(self):
        seq = 'CGCGTGAAAGCAAAACAAGA'
        path = os.path.join(self.base_dir, 'single_seq.fa')
        with btllib.SeqReader(path, btllib.SeqReaderFlag.SHORT_MODE) as reader:
            read = [record.seq for record in reader]
        self.assertListEqual([seq], read)
