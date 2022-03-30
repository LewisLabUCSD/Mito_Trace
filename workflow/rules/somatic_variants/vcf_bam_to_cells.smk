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
allCellBarcodes = {}
for exper in samples_csv:
    curr = pd.read_csv(join(config["params_dir"],samples_csv[exper]))
    for ind, val in curr.iterrows(): #["sample_name"].values:
        sample_name = val["sample_name"]
        name = f"{exper}_{sample_name}"
        allSamples[name] = val["bam_f"]
        allCellBarcodes[name] = val["barcode_f"]


print("allSamples")
print(allSamples)

if "samples" in config:
    allSamples = {x:x for x in config["samples"]}


peak_names = [] if "peaks" not in config else config["peaks"].keys()
peak_names.append("all")

rule all:
    input:
        # expand("{outdir}/aggregate/needle_post/{s}/scPileupVars",
        #     outdir=join(res), s=allSamples.keys()),
        expand("{outdir}/aggregate/needle_post/peaks_{peaks}/{s}/scPileupVars",
            outdir=join(res), s=allSamples.keys(),
               peaks=config.get("peaks", [])+["all"]),



#############################################
# Convert needlestack results to cell-by-vars
############################################
def get_vcf(wildcards):
    if wildcards.peaks == "all":
        return f"{wildcards.outdir}/aggregate/variants.all.vcf"

    return f"{wildcards.outdir}/aggregate/{wildcards.peaks}.variants.vcf"

rule intersect_vcf:
    input:
        bam = "{outdir}/aggregate/preproc/{s}.bam",
        vcf = get_vcf #"{outdir}/aggregate/variants.all.vcf" #vcf = "{outdir}/aggregate/{peaks}.variants.vcf"
    output:
        bam = temp("{outdir}/aggregate/needle_post/peaks_{peaks}/{s}.bam"),
        bai = temp("{outdir}/aggregate/needle_post/peaks_{peaks}/{s}.bam.bai"),
    resources:
        mem_mb=20000
    shell:
        "bedtools intersect -wa -a {input.bam} -b {input.vcf}  > {output.bam} && samtools index {output.bam}"


rule filter_barcodes:
    input:
        bam = "{outdir}/aggregate/needle_post/peaks_{peaks}/{s}.bam",
        bai = "{outdir}/aggregate/needle_post/peaks_{peaks}/{s}.bam.bai",
    output:
        bam = temp("{outdir}/aggregate/needle_post/peaks_{peaks}/cellFilt_{s}.bam"),
        bai = temp("{outdir}/aggregate/needle_post/peaks_{peaks}/cellFilt_{s}.bam.bai"),


def get_barcode(wildcards):
    return allCellBarcodes[wildcards.s]

rule filter_barcodes_01:
    input:
        bam = "{outdir}/aggregate/needle_post/peaks_{peaks}/cellFilt_{s}.bam",
        cells= get_barcode  #"{output}/data/{sample}/MT/cellr_True/{sample}_barcode_data.txt"
        #"{output}/data/{sample}/MT/{sample}.MT.bam",
    output:
        sam = temp("{outdir}/aggregate/needle_post/peaks_{peaks}/filtered.{s}.sam"),
    params:
        cmd = lambda wildcards, input: f"{{ samtools view -H {input.bam} & samtools view {input.bam} | LC_ALL=C grep -F -f {input.cells}; }}"
    resources:
        mem_mb=60000
    run: shell("{params.cmd} > {output.sam} 2> {output.sam}.log")
    #shell:  "{{ samtools view -H {input[0]} & samtools view {input[0]} | LC_ALL=C grep -F -f {input.cells}; }} > {output.sam} 2> {output.sam}.log"


rule filter_barcodes_02:
    input:
        sam = "{outdir}/aggregate/needle_post/peaks_{peaks}/filtered.{s}.sam"
    output:
        bam = "{outdir}/aggregate/needle_post/peaks_{peaks}/filtered.{s}.bam" #temp("{output}/data/{sample}/MT/filtered.bam"),
    shell: "samtools view -b {input.sam} > {output.bam}"

rule sortCB:
    input:
        bam= "{outdir}/aggregate/needle_post/peaks_{peaks}/filtered.{s}.bam",
    output: "{outdir}/aggregate/needle_post/peaks_{peaks}/{s}.CB.bam"  #temp("{output}/data/{sample}/MT/{sample}.MT.CB.bam")
    shell: "samtools sort -t CB {input.bam} > {output} && samtools index {output}"


rule split_by_cb:
    input:
        #bam = "{outdir}/aggregate/needle_post/{s}.bam",
        cb_bam = "{outdir}/aggregate/needle_post/peaks_{peaks}/{s}.CB.bam"
    output:
        scBam = ancient(directory("{outdir}/aggregate/needle_post/peaks_{peaks}/{s}/scBam")),
        #cb_bam = "{outdir}/aggregate/needle_post/peaks_{peaks}/{s}.CB.bam"
    params:
        #script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/needlestack_convert_vars_to_cells.ipynb")
        script = join(ROOT_DIR, "src/mtpreproc/split_by_CB.py")
    resources:
        mem_mb=20000
    shell:
        "python {params.script} {input.bam} {output.scBam}"


rule cb_to_pileup:
    input:
        "{outdir}/aggregate/needle_post/peaks_{peaks}/{s}/scBam"
    output:
        directory("{outdir}/aggregate/needle_post/peaks_{peaks}/{s}/scPileupVars"),
    params:
        script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/scPileup_counts.py"),
        nproc = lambda wildcards: config["nproc"] if "nproc" in config else 24
    threads: 16
    shell: "python {params.script} {input} {output} --nproc {params.nproc}"


rule pileup_to_cell_vars:
    input:
        "{outdir}/aggregate/needle_post/peaks_{peaks}/{s}/scPileupVars",
         get_vcf
        #"{outdir}/aggregate/needle_post.variants.vcf"
    output:
        "{outdir}/aggregate/needle_post/peaks_{peaks}/{s}/cell_by_var/cell_vars.AF.tsv",
        "{outdir}/aggregate/needle_post/peaks_{peaks}/{s}/cell_by_var/cell_vars.DP.tsv"
    params:
        outdir = lambda wildcards, output: dirname(output[0]),
        script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/pileup_to_cell_vars.py")
    shell: "python {params.script} {input}"

