import os
import string
import random
import btllib
import unittest

class SeqWriterTests(unittest.TestCase):

    def test_fasta(self):
        ids = [ "1", "2" ]
        comments = [ "comment1", "comment2" ]
        seqs = [ "ACTG", "TGCA" ]
        quals = [ "!@^&", "(#&$" ]
            
        letters = string.ascii_lowercase
        for iteration in range(3):
            random_filename = ''.join(random.choice(letters) for i in range(64))

            writer_fasta = btllib.SeqWriter(random_filename, btllib.SeqWriter.FASTA)
            for i in range(2):
                writer_fasta.write(ids[i], comments[i], seqs[i], "")
            writer_fasta.close()

            reader_fasta = btllib.SeqReader(random_filename, btllib.SeqReaderFlag.SHORT_MODE)
            self.assertEqual(reader_fasta.get_format(), btllib.SeqReader.SeqReaderFormat_FASTA)

            i = 0
            for record in reader_fasta:
                self.assertEqual(record.id, ids[i])
                self.assertEqual(record.comment, comments[i])
                self.assertEqual(record.seq, seqs[i])
                self.assertTrue(not record.qual)

                i += 1
            self.assertEqual(i, 2)

            reader_fasta.close()
            os.remove(random_filename)
            
    def test_fasta(self):
        ids = [ "1", "2" ]
        comments = [ "comment1", "comment2" ]
        seqs = [ "ACTG", "TGCA" ]
        quals = [ "!@^&", "(#&$" ]
            
        letters = string.ascii_lowercase
        for iteration in range(3):
            random_filename = ''.join(random.choice(letters) for i in range(64))

            writer_fasta = btllib.SeqWriter(random_filename, btllib.SeqWriter.FASTQ)
            for i in range(2):
                writer_fasta.write(ids[i], comments[i], seqs[i], quals[i])
            writer_fasta.close()

            reader_fasta = btllib.SeqReader(random_filename, btllib.SeqReaderFlag.SHORT_MODE)
            self.assertEqual(reader_fasta.get_format(), btllib.SeqReader.SeqReaderFormat_FASTQ)

            i = 0
            for record in reader_fasta:
                self.assertEqual(record.id, ids[i])
                self.assertEqual(record.comment, comments[i])
                self.assertEqual(record.seq, seqs[i])
                self.assertEqual(record.qual, quals[i])

                i += 1
            self.assertEqual(i, 2)

            reader_fasta.close()
            os.remove(random_filename)