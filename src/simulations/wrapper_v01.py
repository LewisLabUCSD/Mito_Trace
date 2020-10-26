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

os.chdir(RESULTS)

default_params_f = join(ROOT_DIR, 'parameters/simulations/simple_v01.yaml')
sweep_params_f = join(ROOT_DIR, 'parameters/simulations/sweep_v01.yaml')
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
