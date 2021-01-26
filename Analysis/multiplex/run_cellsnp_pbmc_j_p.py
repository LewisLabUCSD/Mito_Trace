# Creates cellsnp output files to be used for downstream multiplex.

import snakemake
from src.config import ROOT_DIR
import os

snakefile = os.path.join(ROOT_DIR, "workflows", "deconvolve_samples.Snakefile")
out = snakemake.snakemake(snakefile=snakefile, printshellcmds=True, dryrun=False,
                          verbose=True, configfiles=["parameters/2020_11_18.yaml"],
                          cores=20, resources={"mem_mb":100})
print(out)