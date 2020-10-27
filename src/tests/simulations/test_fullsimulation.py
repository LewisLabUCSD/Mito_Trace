from src.simulations import Simulation, ParameterSweep, FullSimulation
from src.simulations.utils.config import read_config_file
from src.config import RESULTS, ROOT_DIR
import time
from os.path import join
import os
import unittest


os.chdir(RESULTS)
params = os.path.join(ROOT_DIR, 'parameters/simulations/test_simple.yaml')
#params = read_config_file(params)

class TestSum(unittest.TestCase):
    def test_iter(self):
        t = time.time()
        sim = FullSimulation(params)
        sim.run()
        print("time", time.time()-t)
        sim.save()
        return sim

if __name__ == '__main__':
    unittest.main()

