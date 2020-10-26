from src.simulations import Simulation, ParameterSweep, FullSimulation
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
    def test_sweepparams(self):
        default_params_f = join(ROOT_DIR,
                                'parameters/simulations/simple.yaml')
        sweep_params_f = join(ROOT_DIR,
                              'parameters/simulations/test_sweep.yaml')
        sweep = ParameterSweep(default_params_f, sweep_params_f)
        if hasattr(sweep, 'f_save'):
            if os.path.exists(sweep.f_save):
                print('Already ran. Loading')
                sweep.load()
            else:
                sweep.run_sweep()  # Runs the simulation
                sweep.save()
        else:
            print('no attr')
            sweep.run_sweep()  # Runs the simulation
            sweep.save()
        x = sweep.plot_before_after_all_clones()
        self.assertTrue(x is None)

        x = sweep.plot_sensitivity_and_dropout()  # plots results
        self.assertTrue(x is None)


if __name__ == '__main__':
    unittest.main()


