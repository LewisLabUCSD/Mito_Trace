from src.simulations.pipeline import Simulation, FullSimulation
from src.simulations.utils.config import read_config_file
from src.config import RESULTS, ROOT_DIR
import time
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


def testIter():
    t = time.time()
    sim = FullSimulation(params)
    sim.run()
    print("time", time.time()-t)
    sim.save()

    return sim



sim = testIter()

