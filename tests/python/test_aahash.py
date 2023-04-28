import btllib
import unittest

class AAHashTests(unittest.TestCase):
   
    def test_aahash(self):
        seq = "RITMLYTI"
        k = 8
        h = 3

        aahash = btllib.AAHash(seq, h, k, 1)
        
        lvl1_hashes = [
            4971716992320218996, 2124324262552088491, 8599050776293572050,
            4162540294170059916, 4246229847304014199, 3577776033547863068,
            7629354742360239848, 337326900406420216,  4971716992320218996
        ]
        
        num_kmers = len(seq) - k + 2
        idx = 0

        while aahash.roll():
            print(aahash.hashes())
            self.assertEqual(lvl1_hashes[idx], aahash.hashes()[0])
            idx += 1
