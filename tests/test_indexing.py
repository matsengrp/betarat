import unittest
import math
from beta_rat import determine_eval_set, navigate, EvalSet

class TestIndexing(unittest.TestCase):
    def setUp(self):
        self.eval_sets = [EvalSet(math.floor, 0, 0.5) for i in xrange(5)]
        self.navigator = navigate(self.eval_sets, 1000)

    def test_determine_eval_set(self):
        self.assertEqual(determine_eval_set(0, 3), 0)
        self.assertEqual(determine_eval_set(1, 3), 3)
        self.assertEqual(determine_eval_set(2, 3), 2)

    def test_type_of_determine_eval_set(self):
        self.assertIsInstance(determine_eval_set(1, 3), int)

    def test_navigate_object(self):
        nav_order = [self.navigator.next()[1] for i in xrange(16)]
        wanted_order = [self.eval_sets[i] for i in [4,3,4,2,4,3,4,1,4,3,4,2,4,3,4,0]]
        self.assertEqual(nav_order, wanted_order)

    def test_navigate_inner_index(self):
        nav_inner_is = [self.navigator.next()[2] for i in xrange(16)]
        wanted_inner_is = [0,0,1,0,2,1,3,0,4,2,5,1,6,3,7,0]
        self.assertEqual(nav_inner_is, wanted_inner_is)

