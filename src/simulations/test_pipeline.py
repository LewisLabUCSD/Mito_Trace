from src.simulations.pipeline import Simulation
from src.simulations.utils.config import read_config_file
from src.config import RESULTS, ROOT_DIR

import os

os.chdir(RESULTS)
params = os.path.join(ROOT_DIR, 'parameters/simulation/simple.yaml')
params = read_config_file(params)

s = Simulation(params)
s.initialize()
s.grow()

s.save()

# p = s.params
# s.init_clone_dict()
# s.init_cell_coverage()
# s.init_cell_af()