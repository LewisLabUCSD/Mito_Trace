# Creates cellsnp output files to be used for downstream multiplex.

import snakemake
from src.config import ROOT_DIR
import os

snakefile = os.path.join(ROOT_DIR, "workflows", "deconvolve_samples.Snakefile")
out = snakemake.snakemake(snakefile=snakefile, printshellcmds=True,dryrun=False,
                          verbose=True, configfiles=["/data2/mito_lineage/parameters/CHIP_april08_2021/mttrace.yaml"],
                          cores=20, resources={"mem_mb":100}, forcerun=["chrM_enrichment"])
print(out)
