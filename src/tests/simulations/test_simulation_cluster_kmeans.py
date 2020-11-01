from src.simulations import Simulation
from src.simulations.analysis import Analysis
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
    def test_oneIter(self):
        s = Simulation(params)
        s.initialize()
        s.grow()
        s.subsample_new(to_delete=False)

        # Cluster Metrics
        an = Analysis()
        an.cluster_kmeans(s.subsample_new_cell_af)
        # an.cluster_performance()
        # an.estimate_growth_rate()
        # an.compare_befaft_kl()
        # an.compare_kl_across_samples()


        return

if __name__ == '__main__':
    unittest.main()

