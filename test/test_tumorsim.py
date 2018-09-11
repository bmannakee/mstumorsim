# Run tests for tumorsim

import unittest

import mstumorsim.utils as utils
from mstumorsim.tumorsim import SNVtree


class TestTumorSim(unittest.TestCase):
    def test_basic_sim(self):
        tree = SNVtree(1000,empty=False,sigs=[1,5,4,6],timepoints=[0,0,.25,.75])
        tree.run()
        print(f'Number of cells is {len(tree.get_cells())}')
        self.assertGreaterEqual(len(tree.get_cells()),1000)
        self.assertLessEqual(len(tree.get_cells()),1002)

    def test_empty_tree(self):
        tree = SNVtree(1000,empty=True)
        tree.run()
        self.assertGreaterEqual(len(tree.get_cells()),1000)
        self.assertLessEqual(len(tree.get_cells()),1002)

    def test_spectrum(self):
        spec = utils.get_spectrum([1,5])
        print(f'Length of spectrum is {len(spec)}')
        print(f'Sum of spectrum is {sum(spec)}')
        self.assertEqual(len(spec),96)
        self.assertAlmostEqual(sum(spec),1)

    def test_get_mutations(self):
        tree = SNVtree(1000,empty=False,sigs=[1,5,4,6],timepoints=[0,0,.25,.75])
        tree.run()
        m = tree.get_mutations()
        self.assertIsInstance(m,list) 
        # Add code to check tuples once they are written

if __name__ == '__main__':
    unittest.main()