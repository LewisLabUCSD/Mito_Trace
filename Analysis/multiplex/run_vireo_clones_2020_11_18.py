import papermill as pm
from src.config import ROOT_DIR
from src.utils.parse_config import read_config_file
import os


config = read_config_file('parameters/2020_11_18.yaml')

inp = os.path.join(ROOT_DIR, "src", "vireo", "vireoSNP_clones.ipynb" )
if not os.path.exists("data"):
   os.mkdir("data")
outp = "data/vireo_clones/2020_11_18.ipynb"


pm.execute_notebook(
   inp,
   outp,
   parameters = dict(AD_F=config["AD"], DP_F=config["DP"])
)


# Run for MT
#data/{name}_cellSNP_chrM
inp = os.path.join(ROOT_DIR, "src", "vireo", "vireoSNP_donors.ipynb" )
if not os.path.exists("data/vireo_donors"):
   os.mkdir("data/vireo_donors")
outp = "data/vireo_donors/2020_11_18.ipynb"
