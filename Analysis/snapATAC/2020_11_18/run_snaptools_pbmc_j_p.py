import snakemake
from src.config import ROOT_DIR
import os

snakefile = os.path.join(ROOT_DIR, "workflows", "tenx_to_snap.Snakefile")
out = snakemake.snakemake(snakefile=snakefile, verbose=True, configfiles=["snap_pbmc.yaml"], cores=32, resources={"mem_mb":100})
print(out)