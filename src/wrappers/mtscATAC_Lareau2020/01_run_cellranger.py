import snakemake
import os
from src.config import ROOT_DIR, RESULTS

os.chdir(ROOT_DIR)
snakefile = "workflows/cellranger_atac.snakefile"
configfile = "parameters/mtscATAC/Lareau_2020/Lareau_2020.yaml"
#samples = "parameters/mtscATAC/Lareau_2020.csv"

cmd = snakemake.snakemake(snakefile, configfiles=[configfile],
                          report="cellrATAC_lareau.html",
                          verbose=True)
