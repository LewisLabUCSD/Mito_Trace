from os.path import join, dirname, basename
from src.config import ROOT_DIR
from src.utils.parse_config import read_config_file
import os
import numpy as np
from os.path import join, dirname
import pandas as pd
from snakemake.utils import min_version
from icecream import ic
min_version("6.0")
print('config', config)
from src.config import ROOT_DIR

########################################################################
# Setup parameters and outdir
########################################################################
res = join(config["outdir"], config["prefix"])

samples_csv = config["samples_csv"]
allSamples = {}
for exper in samples_csv:
    curr = pd.read_csv(join(config["params_dir"],samples_csv[exper]))
    for ind, val in curr.iterrows(): #["sample_name"].values:
        sample_name = val["sample_name"]
        name = f"{exper}_{sample_name}"
        allSamples[name] = val["bam_f"]

print("allSamples")
print(allSamples)

if "samples" in config:
    allSamples = {x:x for x in config["samples"]}

rule all:
    input:
        expand("{outdir}/aggregate/needle_post/{s}/scPileupVars",
            outdir=join(res), s=allSamples.keys()),

#############################################
# Convert needlestack results to cell-by-vars
#############################################
rule intersect_vcf:
    input:
        bam = "{outdir}/aggregate/preproc/{s}.bam",
        vcf = "{outdir}/aggregate/variants.all.vcf"
    output:
        bam = "{outdir}/aggregate/needle_post/{s}.bam",
        bai = "{outdir}/aggregate/needle_post/{s}.bam.bai",
    resources:
        mem_mb=20000
    shell:
        "bedtools intersect -wa -a {input.bam} -b {input.vcf}  > {output.bam} && samtools index {output.bam}"


rule split_by_cb:
    input:
        bam = "{outdir}/aggregate/needle_post/{s}.bam",
    output:
        scBam = directory("{outdir}/aggregate/needle_post/{s}/scBam"),
        cb_bam = "{outdir}/aggregate/needle_post/{s}.CB.bam"
    params:
        #script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/needlestack_convert_vars_to_cells.ipynb")
        script = join(ROOT_DIR, "src/mtpreproc/split_by_CB.py")
    resources:
        mem_mb=20000
    shell:
        "python {params.script} {input.bam} {output.scBam}"

rule cb_to_pileup:
    input:
        "{outdir}/aggregate/needle_post/{s}/scBam"
    output:
        directory("{outdir}/aggregate/needle_post/{s}/scPileupVars"),
    params:
        script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/scPileup_counts.py"),
        nproc = lambda wildcards: config["nproc"] if "nproc" in config else 24
    threads: 16
    shell: "python {params.script} {input} {output} --nproc {params.nproc}"

rule pileup_to_cell_vars:
    input:
        "{outdir}/aggregate/needle_post/{s}/scPileupVars",
        "{outdir}/aggregate/needle_post.variants.vcf"
    output:
        "{outdir}/aggregate/needle_post/{s}/cell_by_var/cell_vars.AD.tsv",
        "{outdir}/aggregate/needle_post/{s}/cell_by_var/cell_vars.DP.tsv"
    params:
        outdir = lambda wildcards, output: dirname(output[0]),
        script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/pileup_to_cell_vars.py")
    shell: "python {params.script} {input}"

