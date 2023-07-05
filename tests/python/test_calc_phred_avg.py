import btllib
import unittest


class CalcPhredAvgTests(unittest.TestCase):

    def test_calc_phred_avg(self):
        qual = "$$%%)*0)'%%&$$%&$&'''*)(((((()55561--.12356577-++**++,////.*))((()+))**010/..--+**++*+++)++++78883"
        self.assertAlmostEqual(
            btllib.calc_phred_avg(qual, 0, 10), 6.4, places=3)
        self.assertAlmostEqual(btllib.calc_phred_avg(qual), 10.949, places=3)
        self.assertAlmostEqual(
            btllib.calc_phred_avg(qual, 0, 4), 3.5, places=3)
        self.assertAlmostEqual(
            btllib.calc_phred_avg(qual, 5, 20), 6.15, places=3)
