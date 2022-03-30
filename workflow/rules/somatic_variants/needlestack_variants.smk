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

#
rule all:
  input:
    #expand("{outdir}/all_variants.vcf", outdir=join(res, "somatic_variants")),
    expand("{outdir}/aggregate/variants.all.vcf",outdir=res),
    expand("{outdir}/aggregate/{peaks}.variants.vcf",outdir=join(res), peaks=config["bed_regions"].keys()),
    expand("{outdir}/aggregate/preproc/{s}.bam",outdir=join(res), s=allSamples.keys()),
    expand("{outdir}/aggregate/needle_post/{s}/scPileupVars",
        outdir=join(res), s=allSamples.keys()),

def get_bam(wildcards):
    #print(allSamples[wildcards.s])
    return allSamples[wildcards.s]

rule move_bam:
    """ Loop through each run's samples csv files and add to outdir 
    """
    input:
        bam = get_bam
    output:
        bam = "{outdir}/aggregate/preproc/{s}.bam",
        bai = "{outdir}/aggregate/preproc/{s}.bam.bai",
    shell: "ln -s {input.bam} {output.bam} && cp {input.bam}.bai {output.bai}"


# peak files from the merged peaks in annotation
rule needlestack_all:
    input:
        bam_in = expand("{{outdir}}/aggregate/preproc/{s}.bam", s=allSamples.keys()),
        bai_in = expand("{{outdir}}/aggregate/preproc/{s}.bam.bai", s=allSamples.keys()),
        #bed_in = config["merged_peaks_needle_f"]  # "{outdir}/bed/peaks.bed"
    output:
        vcf = "{outdir}/aggregate/variants.all.vcf"
    params:
        bam_indir = lambda wildcards, input: dirname(input.bam_in[0]),
        vcf = lambda wildcards, output: basename(output.vcf),
        ref = join(ROOT_DIR, config["ref_fa"]),
        min_dp = 5,
        nsplit = config["ncpus"]
    shell:
        "cd {wildcards.outdir}/aggregate/ && nextflow run iarcbioinfo/needlestack --input_bams preproc/ --ref {params.ref} --min_dp {params.min_dp} --nsplit {params.nsplit} --use_file_name --output_vcf {params.vcf}"


rule needlestack:
    input:
        bam_in = expand("{{outdir}}/aggregate/preproc/{s}.bam", s=allSamples.keys()),
        bai_in = expand("{{outdir}}/aggregate/preproc/{s}.bam.bai", s=allSamples.keys()),
        #bed_in = config["merged_peaks_needle_f"]  # "{outdir}/bed/peaks.bed"
    output:
        vcf = "{outdir}/aggregate/{peaks}.variants.vcf"
    params:
        bam_indir = lambda wildcards, input: dirname(input.bam_in[0]),
        vcf = lambda wildcards, output: basename(output.vcf),
        ref = join(ROOT_DIR, config["ref_fa"]),
        min_dp = 5,
        bed = lambda wildcards: config["bed_regions"][wildcards.peaks],
        nsplit = config["ncpus"]
    shell:
        "cd {wildcards.outdir}/aggregate/ && nextflow run iarcbioinfo/needlestack --input_bams preproc/ --bed {params.bed} --ref {params.ref} --min_dp {params.min_dp} --nsplit {params.nsplit} --use_file_name --output_vcf {params.vcf}"

    #nextflow run iarcbioinfo/needlestack --input_bams preproc/ --bed /data/Mito_Trace/output/geneRegions/gff_A2/chip_genes.bed --ref /mnt/md0/isshamie/Projects/Mito_Trace/data/processed/genomes/mtMasked/GRCh38_MT_blacklist_A2_2020/fasta/genome.fa --min_dp 15 --nsplit 16 --use_file_name --output_vcf chip_genes.variants.vcf

#
#
# #############################################
# # Convert needlestack results to cell-by-vars
# #############################################
# rule intersect_vcf:
#     input:
#         bam = "{outdir}/aggregate/preproc/{s}.bam",
#         vcf = "{outdir}/aggregate/variants.all.vcf"
#     output:
#         bam = "{outdir}/aggregate/needle_post/{s}.bam",
#         bai = "{outdir}/aggregate/needle_post/{s}.bam.bai",
#     resources:
#         mem_mb=20000
#     shell:
#         "bedtools intersect -wa -a {input.bam} -b {input.vcf}  > {output.bam} && samtools index {output.bam}"
#
#
# rule split_by_cb:
#     input:
#         bam = "{outdir}/aggregate/needle_post/{s}.bam",
#         cb_bam = temp("{outdir}/aggregate/needle_post/{s}.CB.bam"),
#     output:
#         directory("{outdir}/aggregate/needle_post/{s}/scBam")
#     params:
#         #script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/needlestack_convert_vars_to_cells.ipynb")
#         script = join(ROOT_DIR, "src/mtpreproc/split_by_CB.py")
#     resources:
#         mem_mb=20000
#     shell:
#         "python {params.script} {input.bam} {output}"
#
# rule cb_to_pileup:
#     input:
#         "{outdir}/aggregate/needle_post/{s}/scBam"
#     output:
#         directory("{outdir}/aggregate/needle_post/{s}/scPileupVars"),
#     params:
#         script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/scPileup_counts.py"),
#         nproc = lambda wildcards: config["nproc"] if "nproc" in config else 24
#     threads: 16
#     shell: "python {params.script} {input} {output} --nproc {params.nproc}"
#
#
# rule pileup_to_cell_vars:
#     input:
#         "{outdir}/aggregate/needle_post/{s}/scPileupVars",
#         "{outdir}/aggregate/needle_post.variants.vcf"
#     output:
#         "{outdir}/aggregate/needle_post/{s}/cell_by_var/cell_vars.AD.tsv",
#         "{outdir}/aggregate/needle_post/{s}/cell_by_var/cell_vars.DP.tsv"
#     params:
#         outdir = lambda wildcards, output: dirname(output[0]),
#         script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/pileup_to_cell_vars.py")
#     shell: "python {params.script} {input}"
#
# #############################################
# #############################################
#
