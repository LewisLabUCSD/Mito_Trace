from src.simulations import Simulation
from src.simulations.utils.config import read_config_file
from src.config import RESULTS, ROOT_DIR
import time
from os.path import join
import os
import unittest


os.chdir(RESULTS)
params = os.path.join(ROOT_DIR, 'parameters/simulations/simple.yaml')
#params = read_config_file(params)

class TestSum(unittest.TestCase):
    def test_oneIter(self):
        s = Simulation(params)
        s.initialize()
        s.grow()
        s.subsample_new(to_delete=False)
        s.save(f_save='simulation/simple_oneIter.p')
        return

if __name__ == '__main__':
    unittest.main()

