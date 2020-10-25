from src.simulations.pipeline import Simulation, FullSimulation, ParameterSweep
from src.simulations.utils.config import read_config_file
from src.config import RESULTS, ROOT_DIR
import time
from os.path import join
import os

os.chdir(RESULTS)
params = os.path.join(ROOT_DIR, 'parameters/simulation/simple.yaml')
#params = read_config_file(params)


def test_oneIter():
    s = Simulation(params)
    s.initialize()
    s.grow()
    s.subsample_new(to_delete=False)
    s.save(f_save='simulation/simple_oneIter.p')
    return


def test_iter():
    t = time.time()
    sim = FullSimulation(params)
    sim.run()
    print("time", time.time()-t)
    sim.save()
    return sim


def test_sweepparams():
    default_params_f = join(ROOT_DIR, 'parameters/simulation/simple.yaml')
    sweep_params_f = join(ROOT_DIR, 'parameters/simulation/test_sweep.yaml')
    sweep = ParameterSweep(default_params_f, sweep_params_f)
    sweep.run_sweep() # Runs the simulation
    sweep.save()
    sweep.plot_sensitivity_and_dropout() #plots results


test_sweepparams()

