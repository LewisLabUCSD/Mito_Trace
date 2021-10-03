import pathlib
import os
import matplotlib as mpl
mpl.use('agg')

r_kernel = "/usr/bin/Rscript"

############
## Get the project directory
############
path = os.path.abspath(__file__)
dir_path = pathlib.Path(os.path.dirname(path))
ROOT_DIR = dir_path.parents[0]
#print(f"Project Directory: {ROOT_DIR}")


############
# Output directories to be used
############
DATA_DIR = os.path.join(ROOT_DIR, "data")
RESULTS = os.path.join(DATA_DIR, "processed")
#RESULTS_DIR = os.path.join(ROOT_DIR, "models")
FIGURES_DIR = os.path.join(ROOT_DIR, "figures")
NUM_CORES = 32


############
## Call the Pipeline run functions
############
## Example:
# from src import multiplex as ml
# from src import variant as vr
# PIPELINES = {
#     "mttrace": [
#         ("multiplex", ml.run, None),
#         ("msa", vr.run, None),
#         ("couplings", cp.run, None),
#         ("compare", cm.run, None),
#     ]
# }