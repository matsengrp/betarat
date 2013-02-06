import unittest
from beta_rat import h2f1
from scipy import special


class TestH2F1(unittest.TestCase):
    """ Just want to make sure that where the scipy solution is actually computable, it aggrees with our exact
    solution; as well as that where the the former is not computable, that we definitely get something
    different."""
    def test_agreement1(self):
        args = [3, -5, 10, 1.3]
        self.assertAlmostEqual(h2f1(*args), special.hyp2f1(*args))

    def test_agreement2(self):
        args = [5, -8, 20, 0.9]
        self.assertAlmostEqual(h2f1(*args), special.hyp2f1(*args))

    def test_disagreement(self):
        args = [40, -200, 250, 1.3]
        self.assertNotAlmostEqual(h2f1(*args), special.hyp2f1(*args))

