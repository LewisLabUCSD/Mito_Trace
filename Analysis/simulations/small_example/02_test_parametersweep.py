from src.simulations.parametersweep import ParameterSweep
from src.simulations.utils.config import read_config_file
from src.config import RESULTS, ROOT_DIR
import time
from os.path import join
import os
import unittest



def main():
    default_params_f = 'parameters/test_simple.yaml'
    sweep_params_f = 'parameters/test_sweep.yaml'
    sweep = ParameterSweep(default_params_f, sweep_params_f)
    sweep.load_results()
    sweep.plot_before_after_all_clones()
    sweep.plot_sensitivity_and_dropout()  # plots results
    sweep.plot_before_after_true_growth()
    #sweep.save()

#
# class TestSum(unittest.TestCase):
#     def test_sweepparams(self):
#         default_params_f = join(ROOT_DIR,
#                                 'parameters/simulations/test_simple.yaml')
#         sweep_params_f = join(ROOT_DIR,
#                               'parameters/simulations/test_sweep.yaml')
#         sweep = ParameterSweep(default_params_f, sweep_params_f)
#
#         sweep.run_sweep()  # Runs the simulation
#         #sweep.save()




if __name__ == '__main__':
    main()
    #unittest.main()


