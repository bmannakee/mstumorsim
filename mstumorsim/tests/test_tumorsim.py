# Run tests for tumorsim

import unittest

from mstumorsim.tumorsim import SNVtree

class TestTumorSim(unittest.TestCase):
    def test_basic_sim(self):
        spec = [1/96. for x in range(95)]
        spec.append(1-sum(spec))
        tree = SNVtree(100,spec)
        self.assertEqual(len(tree.get_cells,100))
