from src.config import ROOT_DIR, DATA_DIR
import os

os.chdir(ROOT_DIR)
params_f = "parameters/simulation/clones.yaml"


from .dominant_clone import simulate_growth

