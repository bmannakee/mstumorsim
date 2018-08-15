# Run tests for tumorsim

import unittest

import mstumorsim.utils as utils
from mstumorsim.tumorsim import SNVtree


class TestTumorSim(unittest.TestCase):
    def test_basic_sim(self):
        spec = [1/96. for x in range(95)]
        spec.append(1-sum(spec))
        tree = SNVtree(1000,spec)
        tree.run()
        print(f'Number of cells is {len(tree.get_cells())}')
        self.assertEqual(len(tree.get_cells()),1000)

    def test_spectrum(self):
        spec = utils.get_spectrum([1,5])
        print(f'Length of spectrum is {len(spec)}')
        print(f'Sum of spectrum is {sum(spec)}')
        self.assertEqual(len(spec),96)
        self.assertAlmostEqual(sum(spec),1)

if __name__ == '__main__':
    unittest.main()