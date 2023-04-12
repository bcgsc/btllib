import os
import unittest
import btllib

class TestNtHash(unittest.TestCase):

    def test_kmer_hash_values(self):
        seq = "ACATGCATGCA"
        k = 5
        h = 3

        hashes = [
            [0xf59ecb45f0e22b9c, 0x4969c33ac240c129, 0x688d616f0d7e08c3],
            [0x38cc00f940aebdae, 0xab7e1b110e086fc6, 0x11a1818bcfdd553],
            [0x603a48c5a11c794a, 0xe66016e61816b9c4, 0xc5b13cb146996ffe]
        ]

        nthash = btllib.NtHash(seq, h, k)
        ntblind = btllib.BlindNtHash(seq, h, k)

        for expected_hashes in hashes:
            nthash.roll()
            ntblind.roll(seq[ntblind.get_pos() + 1])
            self.assertEqual(expected_hashes, list(nthash.hashes()))
            self.assertEqual(expected_hashes, list(ntblind.hashes()))

    def test_kmer_rolling(self):
        seq = "AGTCAGTC"
        h = 3
        k = 4

        nthash = btllib.NtHash(seq, h, k)
        hashes = []

        while nthash.roll():
            h_vals = list(nthash.hashes())
            hashes.append(h_vals)

        self.assertEqual(len(hashes), len(seq) - k + 1)
        self.assertEqual(hashes[0], hashes[-1])

    def test_kmer_rolling_vs_ntbase_hash_values(self):
        seq = "ACGTACACTGGACTGAGTCT"
        kmer1 = "ACGTACACTGGACTGAGT"
        kmer2 = "CGTACACTGGACTGAGTC"
        kmer3 = "GTACACTGGACTGAGTCT"

        nthash = btllib.NtHash(seq, 3, len(seq) - 2)
        nthash_vector = [
            btllib.NtHash(kmer1, nthash.get_hash_num(), len(kmer1)),
            btllib.NtHash(kmer2, nthash.get_hash_num(), len(kmer2)),
            btllib.NtHash(kmer3, nthash.get_hash_num(), len(kmer3))
        ]

        i = 0
        while nthash.roll() and nthash_vector[i].roll():
            for j in range(nthash.get_hash_num()):
                self.assertEqual(nthash.hashes()[j], nthash_vector[i].hashes()[j])
            i += 1
        self.assertEqual(i, 3)

    def test_canonical_hashing(self):
        seq_f = "ACGTACACTGGACTGAGTCT"
        seq_r = "AGACTCAGTCCAGTGTACGT"
        h = 3

        nthash_f = btllib.NtHash(seq_f, h, len(seq_f))
        nthash_r = btllib.NtHash(seq_r, h, len(seq_r))

        nthash_f.roll()
        nthash_r.roll()
        self.assertEqual(nthash_f.get_hash_num(), nthash_r.get_hash_num())
        self.assertEqual(list(nthash_f.hashes()), list(nthash_r.hashes()))