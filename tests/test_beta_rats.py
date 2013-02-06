import unittest
from beta_rat import BetaRat

class TestSanity(unittest.TestCase):
    def setUp(self):
        self.br = BetaRat(5,5,10,10)

    def test_cdf(self):
        self.assertAlmostEqual(self.br.cdf(1), 0.5, places=4)

    def test_ppf(self):
        self.assertAlmostEqual(self.br.ppf(0.5), 1.0, places=3)


class TestFlippingShit(unittest.TestCase):
    def setUp(self):
        self.br1 = BetaRat(5, 10, 45, 40)
        self.br2 = BetaRat(10, 5, 40, 45)
        self.br2_unflipped = BetaRat(10, 5, 40, 45, no_inverting=True)

    def test_inverted_attrs(self):
        self.assertEqual(self.br2.inverted.a1, 5)
        self.assertEqual(self.br2.inverted.a2, 10)
        self.assertEqual(self.br2.inverted.b1, 45)
        self.assertEqual(self.br2.inverted.b2, 40)
        self.assertTrue(self.br2.inverted)

    def test_inversion_working_on_ppf(self):
        self.assertAlmostEqual(self.br2.ppf(0.5, h_init=1e-5),
                self.br2_unflipped.ppf(0.5, h_init=1e-5),
                places=3)
        self.assertAlmostEqual(self.br2.ppf(.75, h_init=1e-5),
                self.br2_unflipped.ppf(0.75, h_init=1e-5),
                places=3)


