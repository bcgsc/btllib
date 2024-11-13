import btllib
import unittest


class CalcPhredAvgTests(unittest.TestCase):

    def test_calc_phred_avg(self):
        qual = "$$%%)*0)'%%&$$%&$&'''*)(((((()55561--.12356577-++**++,////.*))((()+))**010/..--+**++*+++)++++78883"
        self.assertAlmostEqual(
            btllib.calc_phred_avg(qual, 0, 10), 5.34264, places=3)
        self.assertAlmostEqual(btllib.calc_phred_avg(qual), 8.54327, places=3)
        self.assertAlmostEqual(
            btllib.calc_phred_avg(qual, 0, 4), 3.47128, places=3)
        self.assertAlmostEqual(
            btllib.calc_phred_avg(qual, 5, 20), 5.48923, places=3)
