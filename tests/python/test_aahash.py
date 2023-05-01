import btllib
import unittest


class AAHashTests(unittest.TestCase):

    def test_level_1_aahash(self):
        seq = "RITMLYTI"
        k = 8
        h = 1
        level = 1

        aahash = btllib.AAHash(seq, h, k, level)

        lvl1_hashes = [
            4971716992320218996,
            2124324262552088491,
            8599050776293572050,
            4162540294170059916,
            4246229847304014199,
            3577776033547863068,
            7629354742360239848,
            337326900406420216,
            4971716992320218996
        ]

        num_kmers = len(seq) - k + 2
        idx = 0

        while aahash.roll():
            self.assertEqual(lvl1_hashes[idx], aahash.hashes()[0])
            idx += 1

    def test_level_2_aahash(self):
        seq = "RITMLYTI"
        k = 8
        h = 1
        level = 2

        aahash = btllib.AAHash(seq, h, k, level)

        lvl2_hashes = [
            16984239531076160731,
            2122538798213336131,
            1181634448219258361,
            18364025942200420499,
            15509211063111791705,
            9581438197916650957,
            16088538885945095129,
            6961954824987972818,
            16984239531076160731
        ]

        num_kmers = len(seq) - k + 2
        idx = 0

        while aahash.roll():
            self.assertEqual(lvl2_hashes[idx], aahash.hashes()[0])
            idx += 1

    def test_level_3_aahash(self):
        seq = "RITMLYTI"
        k = 8
        h = 1
        level = 3

        aahash = btllib.AAHash(seq, h, k, level)

        lvl3_hashes = [
            14204475232333016996,
            7411036844653579932,
            4010182074210619839,
            12162032840376646189,
            12359534317481130205,
            12474504577607761213,
            12397957167802074209,
            11500784329128563089,
            14204475232333016996
        ]

        num_kmers = len(seq) - k + 2
        idx = 0

        while aahash.roll():
            self.assertEqual(lvl3_hashes[idx], aahash.hashes()[0])
            idx += 1

    def test_level_2_equivalence(self):
        seq1 = "CGATNDQVWHP"
        seq2 = "CGASNEKIFHP"
        k = 11
        h = 1
        level = 2

        aahash1 = btllib.AAHash(seq1, h, k, level)
        aahash2 = btllib.AAHash(seq2, h, k, level)
        aahash1.roll()
        aahash2.roll()
        self.assertEqual(aahash1.hashes()[0], aahash2.hashes()[0])

    def test_level_3_equivalance(self):
        seq1 = "CGANQVWHP"
        seq2 = "CGTDKIFHP"
        k = 9
        h = 1
        level = 3

        aahash1 = btllib.AAHash(seq1, h, k, level)
        aahash2 = btllib.AAHash(seq2, h, k, level)
        aahash1.roll()
        aahash2.roll()
        self.assertEqual(aahash1.hashes()[0], aahash2.hashes()[0])
