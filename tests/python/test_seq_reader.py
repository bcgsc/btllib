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
        
    def test_fasta(self):
        ids = ["asdf", "ghjk"]
        seqs = ["ACTG", "TGCA"]

        for iteration in range(3):
            reader = btllib.SeqReader(os.path.join(self.base_dir, "../input.fa.gz.bz2.xz"),
                                      btllib.SeqReaderFlag.SHORT_MODE)
            self.assertEqual(reader.get_format(), btllib.SeqReader.SeqReaderFormat_FASTA)

            i = 0
            for record in reader:
                self.assertEqual(record.id, ids[i])
                self.assertEqual(record.seq, seqs[i])
                self.assertEqual(record.qual, "")
                i += 1
            self.assertEqual(i, 2)
