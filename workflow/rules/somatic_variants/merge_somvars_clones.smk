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
import subprocess as subp

########################################################################
# Setup parameters and outdir
########################################################################
# rule all:
#     input:
#         expand("{outdir}/som_dendro_{dendro_thresh}/som_clones.ipynb")


#############################################
# Convert needlestack results to cell-by-vars
############################################
# def get_vcf(wildcards):
#     print(f"{wildcards.outdir}/regions_{wildcards.peaks}/gatk_mutect/variants.vcf.gz")
#     return f"{wildcards.outdir}/regions_{wildcards.peaks}/gatk_mutect/variants.vcf.gz"


######################################################
## Loader scripts for the somatic variant outputs
######################################################
def get_af(wildcards):
    w = wildcards
    if w.somvar_method == "needle":
        samplename = f"{config['experiment']}_{w.sample}"
        return join(config["outdir"], "somatic_variants", config["somvars_dir"], f"aggregate/needle_post/peaks_{w.regions}/{samplename}/cells_vars/vcfpad_1/af.pileup.tsv")
    elif w.somvar_method == "gatk_mutect":
        #curr_dir = join(config["outdir"], "somatic_variants",  f"regions_{w.regions}/gatk_mutect/post/{w.sample}/cells_vars/vcfpad_1")
        return join(config["results"], "somatic_variants",  f"regions_{w.regions}/gatk_mutect/post/pileup/{w.sample}/cells_vars/vcfpad_1/af.pileup.tsv.gz")
    else:
        raise ValueError("somvar_method: needle or mutect")


def get_ref(wildcards):
    w = wildcards
    if w.somvar_method == "needle":
        samplename = f"{config['experiment']}_{w.sample}"
        return join(config["outdir"], "somatic_variants", config["somvars_dir"], f"aggregate/needle_post/peaks_{w.regions}/{samplename}/cells_vars/vcfpad_1/af.ref.pileup.tsv")
    elif w.somvar_method == "gatk_mutect":
        return join(config["results"], f"somatic_variants/regions_{w.regions}/gatk_mutect/post/pileup/{w.sample}/cells_vars/vcfpad_1/af.ref.pileup.tsv.gz")
    else:
        raise ValueError("somvar_method: needle or mutect")


def vars_anno(wildcards):
    w = wildcards
    if w.somvar_method == "needle":
        return join(config["outdir"], "somatic_variants", config["somvars_dir"], f"aggregate/needle_post/peaks_{w.regions}/variants.annotate.gene.vcf")
    elif w.somvar_method == "gatk_mutect":
        return join(config["results"], f"somatic_variants/regions_{w.regions}/gatk_mutect/post/variants.annotate.gene.vcf")
    else:
        raise ValueError("somvar_method: needle or mutect")
    return


rule merge_som_clones:
    input:
        af = get_af,
        ref = get_ref,
        cells_meta = "{outdir}/cells_meta.tsv",
        vars_f = vars_anno,
        #clones = "{outdir}/barcodes/_clone_complete.txt",
        #ref= "{out_dir}/aggregate/needle_post/peaks_{peaks}/{sample}/cells_vars/vcfpad_1/af.ref.pileup.tsv",
        #af= "{out_dir}/aggregate/needle_post/peaks_{peaks}/{sample}/cells_vars/vcfpad_1/af.pileup.tsv"
    output:
          note="{outdir}/somatic_variants/{somvar_method}/peaks_{regions}/{sample}/som_clones.ipynb"
          # figs=multiext("{outdir}/somatic_variants/needle/peaks_{regions}/{sample}/",
          #          "clone_altAndRef_numCells.png","clone_altAndRef_sumCounts.png", "donor_altAndRef_numCells.png","donor_altAndRef_sumCounts.png")
    params:
        indir = lambda wildcards, input: dirname(input.af),
        outdir = lambda wildcards, output: dirname(output.note),
        #cells_meta = lambda wildcards, input: join(dirname(dirname(input.clones)), "cells_meta.tsv"),
        script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/vars_to_clones.ipynb")
    resources:
        mem_mb=80000
    shell: "papermill -p cells_meta_f {input.cells_meta} -p indir {params.indir} -p condition {wildcards.sample} -p outdir {params.outdir} -p vars_f {input.vars_f} {params.script} {output.note}"

#
# rule som_clones_samples:
#     input:
#          expand("{{outdir}}/somatic_variants/{{somvar_method}}/peaks_{{regions}}/{sample}/som_clones.ipynb",
#              sample = config["samples"].index)
#          #note="{outdir}/somatic_variants/{somvar_method}/peaks_{regions}/{sample}/som_clones.ipynb",
#     output:
#           vars="{outdir}/somatic_variants/{somvar_method}/peaks_allSamples_{regions}/allSamples.variants.tsv",
#           stats="{outdir}/somatic_variants/{somvar_method}/peaks_allSamples_{regions}/allSamples.variants.stats.tsv"
#     run:
#         allVars = []
#         for i in input:
#             curr_dir = dirname(i)
#             curr_f = join(curr_dir, "clones_alt_numCells.csv")
#             if os.path.exists(curr_f):
#                 curr_vars = pd.read_csv(curr_f)
#                 allVars.append(curr_vars)
#         if len(allVars) != 0:
#             pd.concat(allVars).to_csv(output.vars)
#         else:
#             os.system(f"touch {output.vars}")

rule som_dendro:
    input:
        note = expand("{{outdir}}/somatic_variants/{{somvar_method}}/peaks_{{regions}}/{sample}/som_clones.ipynb",
             sample = config['samples'].index),
        vars_f = vars_anno, #"{outdir}/regions_{peaks}/gatk_mutect/post/variants.annotate.gene.vcf"
        dendro_f = expand("{{outdir}}/barcodes/btwnClones_dendro_dt_{{dendro_thresh}}/donor{d}.mean.csv",
                          d= np.arange(config["N_DONORS"])),
         #"{outdir}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.mean.csv",
        cells_meta = "{outdir}/cells_meta.tsv",
    output:
        note="{outdir}/somatic_variants/{somvar_method}/peaks_{regions}/som_dendro_{dendro_thresh}/som_clones.ipynb",
    params:
        indir = lambda wildcards, input: dirname(dirname(input.note[0])),
        outdir = lambda wildcards, output: dirname(output.note),
        dendro_indir = lambda wildcards, input: dirname(input.dendro_f[0]),
        #cells_meta = lambda wildcards, input: join(dirname(dirname(input.clones)), "cells_meta.tsv"),
        script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/dendro_vars_to_clones.ipynb")
    shell: "papermill -p cells_meta_f {input.cells_meta} -p indir {params.indir} -p dendro_indir {params.dendro_indir} -p outdir {params.outdir} -p vars_f {input.vars_f} {params.script} {output.note}"


rule chip_som_dendro:
    input:
        note = expand("{{outdir}}/somatic_variants/{{somvar_method}}/peaks_{{regions}}/{sample}/som_clones.ipynb",
             sample = config['samples'].index),
        vars_f = vars_anno, #"{outdir}/regions_{peaks}/gatk_mutect/post/variants.annotate.gene.vcf"
        dendro_f = expand("{{outdir}}/barcodes/btwnClones_dendro_dt_{{dendro_thresh}}/donor{d}.mean.csv",
                          d= np.arange(config["N_DONORS"])),
         #"{outdir}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.mean.csv",
        cells_meta = "{outdir}/cells_meta.tsv",
    output:
        note="{outdir}/somatic_variants/{somvar_method}/peaks_{regions}/filt_chip/som_dendro_{dendro_thresh}/som_clones.ipynb",
    params:
        indir = lambda wildcards, input: dirname(dirname(input.note[0])),
        chip_genes = config["bed_regions"]["chip_genes"],
        outdir = lambda wildcards, output: dirname(output.note),
        dendro_indir = lambda wildcards, input: dirname(input.dendro_f[0]),
        #cells_meta = lambda wildcards, input: join(dirname(dirname(input.clones)), "cells_meta.tsv"),
        script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/dendro_vars_to_clones.ipynb")
    shell: "papermill -p cells_meta_f {input.cells_meta} -p indir {params.indir} -p dendro_indir {params.dendro_indir} -p chip_genes_f {params.chip_genes}  -p outdir {params.outdir} -p vars_f {input.vars_f} {params.script} {output.note}"




def get_output(wildcards):
    w = wildcards
    return


rule summary_som:
    input:
        note="{outdir}/somatic_variants/{somvar_method}/peaks_{regions}/som_dendro_{dendro_thresh}/som_clones.ipynb",
        cells_meta = "{outdir}/cells_meta.tsv",
        dendro_f = expand("{{outdir}}/barcodes/btwnClones_dendro_dt_{{dendro_thresh}}/donor{d}.mean.csv",
                  d= np.arange(config["N_DONORS"])),
    output:
        note="{outdir}/somatic_variants/{somvar_method}/peaks_{regions}/som_dendro_{dendro_thresh}/clones_variants_summary/summary_som.ipynb",
    params:
        indir = lambda wildcards, input: dirname(dirname(input.note)),
        outdir = lambda wildcards, output: dirname(output.note),
        dendro_indir = lambda wildcards, input: dirname(input.dendro_f[0]),
        script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/summarize_vars_clones.ipynb"),
    shell: "papermill -p cells_meta_f {input.cells_meta} -p indir {params.indir} -p outdir {params.outdir} {params.script} {output.note}"


rule chip_summary_som:
    input:
        note="{outdir}/somatic_variants/{somvar_method}/peaks_{regions}/som_dendro_{dendro_thresh}/som_clones.ipynb",
        cells_meta = "{outdir}/cells_meta.tsv",
        dendro_f = expand("{{outdir}}/barcodes/btwnClones_dendro_dt_{{dendro_thresh}}/donor{d}.mean.csv",
                  d= np.arange(config["N_DONORS"])),
    output:
        note="{outdir}/somatic_variants/{somvar_method}/peaks_{regions}/som_dendro_{dendro_thresh}/chip_clones_variants_summary/summary_som.ipynb",
    params:
        chip_genes = config["bed_regions"]["chip_genes"],
        indir = lambda wildcards, input: dirname(dirname(input.note)),
        outdir = lambda wildcards, output: dirname(output.note),
        dendro_indir = lambda wildcards, input: dirname(input.dendro_f[0]),
        script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/summarize_vars_clones.ipynb"),
    shell: "papermill -p cells_meta_f {input.cells_meta} -p indir {params.indir} -p dendro_indir {params.dendro_indir} -p outdir {params.outdir} -p to_chip True -p chip_genes_f {params.chip_genes} {params.script} {output.note}"

