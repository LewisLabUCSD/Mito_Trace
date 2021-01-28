# Creates cellsnp output files to be used for downstream multiplex.

import snakemake
from src.config import ROOT_DIR
import os

snakefile = os.path.join(ROOT_DIR, "workflows", "deconvolve_samples.Snakefile")
out = snakemake.snakemake(snakefile=snakefile, printshellcmds=True, dryrun=False,
                          verbose=True, configfiles=["parameters/jan21_2021.yaml"],
                          cores=20, resources={"mem_mb":100})
print(out)