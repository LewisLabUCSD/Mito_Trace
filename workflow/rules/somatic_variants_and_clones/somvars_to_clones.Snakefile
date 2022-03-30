import os
import numpy as np
from os.path import join, dirname
import pandas as pd
from snakemake.utils import min_version
from icecream import ic
min_version("6.0")
print('config', config)
from src.config import ROOT_DIR


rule pileup_to_cell_vars:
    input:
        "{outdir}/aggregate/needle_post/{s}/scPileupVars",
        "{outdir}/aggregate/needle_post.variants.vcf"
    output:
        "{outdir}/aggregate/needle_post/{s}/cell_by_var/cell_vars.AF.tsv",
        "{outdir}/aggregate/needle_post/{s}/cell_by_var/cell_vars.DP.tsv"
    params:
        outdir = lambda wildcards, output: dirname(output[0]),
        script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/pileup_to_cell_vars.py")
    shell: "python {params.script} {input}"


rule merge_som_clones:
    input:
        clones = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/_clone_complete.txt",
        AF = "{outdir}/aggregate/needle_post/{s}/cell_by_var/cell_vars.AF.tsv",
        DP = "{outdir}/aggregate/needle_post/{s}/cell_by_var/cell_vars.DP.tsv",
    output:
          multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/",
                   "clones_somvars.png", "somvars_donors.png", "somvars_clones.tsv")

    params:
        script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/pileup_to_cell_vars.py")
