import unittest
from betarat import BetaRat

class TestSanity(unittest.TestCase):
    def setUp(self):
        self.br = BetaRat(5,5,10,10)

    def test_cdf(self):
        self.assertAlmostEqual(self.br.cdf(1), 0.5, places=4)

    def test_ppf(self):
        self.assertAlmostEqual(self.br.ppf(0.5), 1.0, places=3)

class TestHandlingZeros(unittest.TestCase):
    """ Having some issues with the new priors leading to Zero divisions. Gonna try to fix this up."""
    def test_creation(self):
        BetaRat(0, 0, 0, 0, prior=(0,0))

    def test_pdf(self):
        BetaRat(0, 0, 0, 0, prior=(0,0)).pdf(0.5)
        self.assertEqual(BetaRat(0, 0, 0, 0, prior=(0,0)).pdf(0.0), 0)

    def test_pdf_on_0(self):
        self.assertEqual(BetaRat(1, 3, 5, 7).pdf(0.0), 0)

class TestFlippingShit(unittest.TestCase):
    # Only makes sense for simpson...
    def setUp(self):
        self.br1 = BetaRat(5, 10, 45, 40, prior=(0,0))
        self.br2 = BetaRat(10, 5, 40, 45, prior=(0,0))
        self.br2_unflipped = BetaRat(10, 5, 40, 45, no_inverting=True, prior=(0,0))

    def test_inverted_attrs(self):
        self.assertEqual(self.br2.inverted.a1, 5)
        self.assertEqual(self.br2.inverted.a2, 10)
        self.assertEqual(self.br2.inverted.b1, 45)
        self.assertEqual(self.br2.inverted.b2, 40)
        self.assertTrue(self.br2.inverted)

    def test_inversion_working_on_ppf(self):
        print "here we are"
        self.assertAlmostEqual(self.br2.ppf(0.5, h_init=1e-5, method="simpson"),
                self.br2_unflipped.ppf(0.5, h_init=1e-5, method="simpson"),
                places=3)
        self.assertAlmostEqual(self.br2.ppf(.75, h_init=1e-5, method="simpson"),
                self.br2_unflipped.ppf(0.75, h_init=1e-5, method="simpson"),
                places=3)


