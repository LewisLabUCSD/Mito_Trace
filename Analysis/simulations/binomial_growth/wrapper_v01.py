"""
Running the first version of the lineage tracing simulation.
Files are simple_v01.yaml and wrapper_v01.yaml.
1000 iterations, 10000 num_cells
"""

from src.simulations import Simulation, ParameterSweep, FullSimulation
from src.simulations.utils.config import read_config_file
from src.config import RESULTS, ROOT_DIR
import time
from os.path import join
import os

#os.chdir(RESULTS)

default_params_f = 'parameters/simple_v01.yaml'
sweep_params_f = 'parameters/sweep_v01.yaml'
sweep = ParameterSweep(default_params_f, sweep_params_f)

sweep.run_sweep()  # Runs the simulation
sweep.save()
sweep.plot_before_after_all_clones()
sweep.plot_sensitivity_and_dropout()  # plots results