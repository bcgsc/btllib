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

    def test_kmer_back_rolling(self):
        seq = "ACTAGCTG"
        h = 3
        k = 5

        nthash = btllib.NtHash(seq, h, k)
        hashes = []

        while nthash.roll():
            h_vals = list(nthash.hashes())
            hashes.append(h_vals)

        self.assertEqual(len(hashes), len(seq) - k + 1)

        while hashes:
            self.assertEqual(list(nthash.hashes()), hashes.pop())
            nthash.roll_back()


    def test_skipping_ns(self):
        seq = "ACGTACACTGGACTGAGTCT"
        seq_with_ns = seq

        self.assertGreaterEqual(len(seq_with_ns), 10)
        middle = len(seq_with_ns) // 2
        seq_with_ns = seq_with_ns[:middle] + "NN" + seq_with_ns[middle + 2:]
        k = (len(seq) - 2) // 2 - 1
        nthash = btllib.NtHash(seq_with_ns, 3, k)

        positions = list(range(0, len(seq_with_ns) // 2 - k + 1)) + \
                    list(range(len(seq_with_ns) // 2 + 2, len(seq_with_ns) - k + 1))

        i = 0
        while nthash.roll():
            self.assertEqual(nthash.get_pos(), positions[i])
            i += 1
        self.assertEqual(len(positions), i)
 
    def test_rna(self):
        seq = "ACGTACACTGGACTGAGTCT"
        dna_nthash = btllib.NtHash(seq, 3, 20)

        rna_seq = "ACGUACACUGGACUGAGUCU"
        rna_nthash = btllib.NtHash(rna_seq, 3, 20)

        dna_nthash.roll()
        rna_nthash.roll()

        for i in range(dna_nthash.get_hash_num()):
            self.assertEqual(dna_nthash.hashes()[i], rna_nthash.hashes()[i])
    
    def test_spaced_seed_hash_values(self):
        seq = "ACATGCATGCA"
        seeds = ["11100111"]
        s_seeds = btllib.parse_seeds(seeds)
        k = len(seeds[0])
        h = 3

        true_hashes = [
            [0x10be4904ad8de5d, 0x3e29e4f4c991628c, 0x3f35c984b13feb20],
            [0x8200a7aa3eaf17c8, 0x344198402f4c2a9c, 0xb6423fe62e69c40c],
            [0x3ce8adcbeaa56532, 0x162e91a4dbedbf11, 0x53173f786a031f45]
        ]

        nthash = btllib.SeedNtHash(seq, s_seeds, h, k)

        for h_vals in true_hashes:
            nthash.roll()
            self.assertEqual(list(nthash.hashes()), h_vals)
            
    def test_spaced_seeds(self):
        seq = "ACGTACACTGGACTGAGTCT"
        seeds = ["111110000000011111", "111111100001111111"]
        s_seeds = btllib.parse_seeds(seeds)

        seqM1 = "ACGTACACTTGACTGAGTCT"
        seqM2 = "ACGTACACTGTACTGAGTCT"
        seqM3 = "ACGTACACTGCACTGAGTCT"

        k = len(seq) - 2
        self.assertEqual(k, len(seeds[0]))
        self.assertEqual(k, len(seeds[1]))

        seed_nthash = btllib.SeedNtHash(seq, s_seeds, 2, k)
        seed_nthashM1 = btllib.SeedNtHash(seqM1, s_seeds, 2, k)
        seed_nthashM2 = btllib.SeedNtHash(seqM2, s_seeds, 2, k)
        seed_nthashM3 = btllib.SeedNtHash(seqM3, s_seeds, 2, k)

        hashes = []

        self.assertEqual(seed_nthash.get_hash_num(), len(seeds) * 2)

        steps = 0
        while seed_nthash.roll():
            self.assertTrue(seed_nthashM1.roll())
            self.assertTrue(seed_nthashM2.roll())
            self.assertTrue(seed_nthashM3.roll())

            seq_sub = seq[steps:steps + k]
            seqM1_sub = seqM1[steps:steps + k]
            seqM2_sub = seqM2[steps:steps + k]
            seqM3_sub = seqM3[steps:steps + k]

            seed_nthash_base = btllib.SeedNtHash(seq_sub, s_seeds, 2, k)
            seed_nthashM1_base = btllib.SeedNtHash(seqM1_sub, s_seeds, 2, k)
            seed_nthashM2_base = btllib.SeedNtHash(seqM2_sub, s_seeds, 2, k)
            seed_nthashM3_base = btllib.SeedNtHash(seqM3_sub, s_seeds, 2, k)

            self.assertTrue(seed_nthash_base.roll())
            self.assertTrue(seed_nthashM1_base.roll())
            self.assertTrue(seed_nthashM2_base.roll())
            self.assertTrue(seed_nthashM3_base.roll())

            hashes.append([])
            for i in range(seed_nthash.get_hash_num()):
                hval = seed_nthash.hashes()[i]
                self.assertEqual(hval, seed_nthashM1.hashes()[i])
                self.assertEqual(hval, seed_nthashM2.hashes()[i])
                self.assertEqual(hval, seed_nthashM3.hashes()[i])
                self.assertEqual(hval, seed_nthashM1_base.hashes()[i])
                self.assertEqual(hval, seed_nthashM2_base.hashes()[i])
                self.assertEqual(hval, seed_nthashM3_base.hashes()[i])
                hashes[-1].append(hval)

            if seed_nthash.get_pos() > 0:
                seed_nthash.peek_back()
                for i in range(seed_nthash.get_hash_num()):
                    self.assertEqual(seed_nthash.hashes()[i], hashes[-2][i])

                seed_nthash.peek_back(seq[seed_nthash.get_pos() - 1])
                for i in range(seed_nthash.get_hash_num()):
                    self.assertEqual(seed_nthash.hashes()[i], hashes[-2][i])

            steps += 1

        self.assertFalse(seed_nthashM1.roll())
        self.assertFalse(seed_nthashM2.roll())
        self.assertFalse(seed_nthashM3.roll())
        self.assertEqual(steps, len(seq) - k + 1)

    def test_base_substitution(self):
        seq = "ACGTACACTGGACTGAGTCT"
        sub = "ACGCGCACTGGACTGAGTCT"

        nthash = btllib.NtHash(seq, 3, len(seq))
        nthash_subbed = btllib.NtHash(sub, 3, len(sub))

        nthash.roll()
        nthash.sub([3, 4], [ord('C'), ord('G')]) ## convert to ASCII
        nthash_subbed.roll()

        self.assertEqual(nthash.get_hash_num(), nthash_subbed.get_hash_num())
        for i in range(nthash.get_hash_num()):
            self.assertEqual(nthash.hashes()[i], nthash_subbed.hashes()[i])

    def test_kmer_pos(self):
        seq = "CCCTATTAGTACAGTAGTGCCTTCATCGGC"
        h = 3
        k = 6
        nthash = btllib.NtHash(seq, h, k)
        
        for index in range(len(seq) - k):
            nthash.roll()
            self.assertEqual(nthash.get_pos(),index)

def test_kmer_peeking(self):
        seq = "ACTGATCAG"
        h = 3
        k = 6

        nthash = btllib.NtHash(seq, h, k)
        nthash.roll()

        steps = 3
        while steps:
            h_current = list(nthash.hashes())
            _ = nthash.peek(seq[nthash.get_pos() + k])
            h_peek_next = list(nthash.hashes())
            self.assertEqual(h_current, h_peek_next)
            nthash.roll()
            self.assertEqual(list(nthash.hashes()), h_current)
            steps -= 1
