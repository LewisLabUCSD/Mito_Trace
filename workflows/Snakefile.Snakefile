#from src.config import ROOT_DIR
#workdir: ROOT_DIR
wildcard_constraints:
    cellr='True|False'

from src.config import ROOT_DIR
from src.utils.parse_config import read_config_file
import os
import numpy as np
from os.path import join, dirname
import pandas as pd
from snakemake.utils import min_version
from icecream import ic
min_version("6.0")


# Setup parameters and outdir
#res = config["results"]
res = join(config["outdir"], "pipeline", config["prefix"])

params = read_config_file(config["config"])
samples = pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)

# Merge the experiment-specific config with the pipeline parameters
params["results"] = res
#params["samples_meta"] = config["samples_meta"]
params["samples"] = samples
params["genome"] = config["genome"]
params["N_DONORS"] = config["N_DONORS"]
####


params_mt = params["mtpreproc"]
params_ft = params["filters"]
ft = params_ft["params"]
params_mult = params["multiplex"]
params_mgatk = params["mgatk"]
params_clones = params["clones"]

cellrbc = params_mt["params"]["cellrbc"]
num_reads_filter = params_mt["params"]["numreadsfilter"]
maxBP = params_mt["maxBP"]
ref_fa = params["genome_path"][config["genome"]]["ref_fa"]

ncellsthreshmgatk = params_mgatk["params"]["ncellsthresh"]
num_cells = params_mult["pseudo_multiplex"]["num_cells"]
is_prop = params_mult["pseudo_multiplex"]["is_proportional"]


# def get_ref(wildcards):
#     w = wildcards
#     mito = config['mito_character']
#     return f"{w.results}/{w.sample}/mapq_{w.mapq}/cellr_{w.cellrbc}/{w.sample}_{w.num_read}/{mito}_refAllele.txt"

#workdir: config["work_dir"]
rule all:
    input:
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/clones/variants_{variants}/{method}/.completed",
        #        output=res, cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], method=params_clones["method"], variants=params_clones["variants"]),
        expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/.pipeline",
            output=res, cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], method=params_clones["method"], variants=params_clones["variants"]),
         # A. Preproc: Initial MT coverage per position per cell
        expand("{output}/figures/{sample}/MT/cellr_{cellrbc}/numread_{num_read}_MT_position_coverage.png",
               output=res,sample=samples.index, num_read=num_reads_filter, cellrbc=cellrbc),
         # B. Multiplex: Multiplexed VAF and depth dendrograms for each donor and selected variants
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/dendrograms/figures/donor{d}_dendrogram.png",
        #        output=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #        mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #        hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], d=np.arange(config["N_DONORS"])),
        expand("{output}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/donor{d}.pdf",
               output=res,cellrbc=cellrbc, num_read=num_reads_filter,
               mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
               hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], d=np.arange(config["N_DONORS"])),
        expand("{output}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/multiplex.pdf",
               output=res,cellrbc=cellrbc, num_read=num_reads_filter,
               mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
               hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], d=np.arange(config["N_DONORS"])),

        # # C. Clones: Variant types
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb",
        #        output=res, cellrbc=cellrbc, num_read=num_reads_filter,
        #        mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #        hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], n_clones=config['multiplex']['n_clone_list']),
        # # D. Enrichment: Volcano plot
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/enrichment/.status",
        #        output=res, cellrbc=cellrbc, num_read=num_reads_filter,
        #        mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #        hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh']),
        # # E. Lineage labels for cells, which can be input to loupe browser
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/clones/variants_{variants}/{method}_{params_name}{method_params}/donor{d}/cells_BC.csv",
        #        d=range(config["N_DONORS"]),
        #        method=config["clones"]["method"],
        #        method_params=[config[x]['params'] for x in config["method"]])
        #

################################################################
## Import from prior snakefile modules
## Here, we redefine the input to be based on our config['files'] dictionary
from snakemake.utils import min_version
min_version("6.0")
module mtpreprocMod:
    snakefile: "./rules/mt_preprocess.smk"
    config: params

module mgatkMod:
    snakefile: "./rules/mgatk.smk"
    config: params

module multMod:
    snakefile: "./rules/multiplex.smk"
    config: params

module lineageMod:
    snakefile: "./rules/lineage.smk"
    config: params




## 1. Go from 10x output to filtered scPileup outdir.
use rule * from mtpreprocMod

# use rule plot_sc_coverageBar_and_heat from mtpreprocMod with:
#     output:
#         save_f_coverage="{output}/figures/{sample}/MT/cellr_{cellrbc}/numread_{num_read}_MT_position_coverage.png",
#         save_f_heat="{output}/figures/{sample}/MT/cellr_{cellrbc}/numread_{num_read}_MT_position.png"


# use rule create_filters from mtpreprocMod as cf with:
#     output:
#         cov = "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/{sample}.coverage.txt",
#         af = "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/af_by_cell.tsv"


## 2. Call variants using MGATK and convert to the Vireo multiplex input
use rule * from mgatkMod
# use rule get_refAllele from mgatkMod with:
#     output:
#         "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/chrM_refAllele.txt"
#
#
# use rule mgatk from mgatkMod with:
#     input:
#         all = rules.create_filters.output,
#         refAllele = rules.get_refAllele.output[0] #"{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/chrM_refAllele.txt"
#     output:
#         "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/{sample}/{sample}.variant.rds"

use rule mgatk_to_vireoIn from mgatkMod with:
    input:
        "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/{sample}.variant.rds"
    output:
        "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/{sample}/vireoIn/cellSNP.tag.AD.mtx"



## 3. Merge the sample vireo inputs together by concatenating.
##### The cell barcodes get changed here!
rule merged:
    """ Merge the sample pileup matrices
    """
    input:
        lambda wildcards: expand("{{output}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/{sample}/vireoIn/cellSNP.tag.AD.mtx",
                                        sample=samples.index)# , minC=wildcards.minC, minAF=wildcards.minAF)
    output: "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/cellSNP.tag.AD.mtx"
    #input: lambda wildcards: expand("data/{{prefix}}/chrM/{name}_cellSNP_minC{{mt_minC}}_minAF{{mt_minAF}}", name=config["samples"])
    params:
        num_cells=num_cells,
        is_prop = is_prop,
        indirs = lambda wildcards, input: [dirname(x) for x in input],
        outdir = lambda wildcards, output: dirname(output[0]),
        prefix= ','.join(samples.index)
    #log: "logs/{prefix}/vireo/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}.log"

    shell: "python -m src.pseudo_batch {params.outdir} {params.indirs} --num_cells {params.num_cells} --is_prop {params.is_prop} --samples {params.prefix}" # > {log} 2>&1"



## 4. De-multiplex the merged conditions using Vireo
use rule * from multMod

use rule donors_type_variants from multMod with:
    output:
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/variants/variants.ipynb",

use rule donors_plotAF from multMod with:
    input:
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb",
    output:
        expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/dendrograms/figures/donor{d}_dendrogram.{suf}",
               d=np.arange(config["N_DONORS"]), suf=["png", "depth.png", "withHigh.png"]),
        #"{output}/merged/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/dendrograms/af_dendro.ipynb",


## 5. Detect clones, depending on the method parameter
#use rule * from lineageMod
nclonelist = params_clones['vireo']['params']['nclonelist']

# use rule vireo_donor_mgatk from lineageMod with:
#     input:
#         mtx=expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/nclones{{nclones}}/donor{d}/cellSNP.tag.AD.mtx",
#                d=np.arange(config["N_DONORS"])),
#         cells=expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/nclones{{nclones}}/donor{d}/cells_meta.tsv",
#                d=np.arange(config["N_DONORS"])),
#     output:
#         report("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/nclones{nclones}/donor{d}.labels.png"),
#         report("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/nclones{nclones}/donor{d}.variants.labels.png"),
#
# use rule * from lineageMod
#
# use rule enrichment_vireo from lineageMod with:
#     input:
#         expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{{variants}}/vireo/nclones{{nclones}}/donor{d}/cells_meta.tsv",
#                d=np.arange(config["N_DONORS"]))
#     output:
#         "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/enrichment/volcano_Fisher_foldNorm.png"
#
#
#
# use rule complete_lineage from lineageMod with:
#     input:
#         "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/temp/.tmp",
#     output:
#         "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/.completed"
#



rule scPileup_filter_mgatk:
    input:
        cov =  expand("{{output}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/{sample}.coverage.txt",
                      sample=samples.index),
        mult = "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb",
    output:
        cov = expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/d{d}.coverage.txt",
                     d=range(config["N_DONORS"]))
    params:
        #outdir = lambda wildcards, output: dirname(output.cov),
        cells_meta = lambda wildcards, input: join(dirname(input.mult), "cells_meta.tsv"),
        sample = samples['sample_name'].values
    script: join(ROOT_DIR, "src/donor_filter_mgatk.py")


# rule donor_get_refAllele:
#     #input: config["mt_ref_fa"],
#     params: params["mgatk"]["chrM_refAllele"]
#     output: "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/chrM_refAllele.txt"
#     shell: 'cp {params} {output}'
#

# use rule mgatk from mgatkMod with:
#     input:
#         all = "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/d{d}.coverage.txt",
#         refAllele = "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/chrM_refAllele.txt"
#     output:
#         vars_f = "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/mgatk_donor/d{d}.variant.rds",
#         vars_qc = "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/mgatk_donor/d{d}.variantQC.png"
#     params:
#         sample = lambda wildcards: f"d{wildcards.d}",
#         outdir =  "donor_mgatk",
#         data_dir =  lambda wildcards, input: dirname(input.all)

rule donor_mgatk:
    input:
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/d{d}.coverage.txt",
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/chrM_refAllele.txt"
    output:
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/mgatk_donor/d{d}.variant.rds",
        report("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/mgatk_donor/d{d}.variantQC.png", category="Variants", subcategory="donorMGATK")
    params:
        data_dir =lambda wildcards, input: dirname(input[0]),
        donor = lambda wildcards: f"d{wildcards.d}",
    shell:
        "./R_scripts/wrap_mgatk.R {params.data_dir} mgatk_donor/{params.donor} FALSE"


########################################################################
## Workflow A: Directly from multiplex output
########################################################################
rule vireo:
    input:
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb" #get_input_multiplex
    output:
        expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{{nclones}}/donor{d}.labels.png", d=np.arange(config["N_DONORS"])),#, category="lineage"),
        expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{{nclones}}/donor{d}.variants.labels.png", d=np.arange(config["N_DONORS"])),#, category="lineage"),
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{nclones}/cells_meta.tsv"
    params:
        notebook=join("src", "vireo", "2_MT_Lineage_Construct.ipynb"),
        #output_notebook = lambda wildcards, output: join(dirname(output[0]), 'clones.ipynb'),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["N_DONORS"],
        nclones= lambda wildcards: wildcards.nclones #",".join([str(x) for x in nclonelist])
    threads: 8
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p nclones {params.nclones} {params.notebook} {params.OUTDIR}/output.ipynb"


########################################################################
## Workflow B: After multiplexing, separate by donor, grab variants from filters
## that overlap with both conditions, and then run mgatk to call variants again.
##
# Extract pileups for each donor from original filter matrices
########################################################################
rule donor_mgatk_to_vireoIn:
    input: "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/mgatk_donor/d{d}.variant.rds"
    output:
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/donor{d}/mgatk_donor/cellSNP.tag.AD.mtx",
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        donor = lambda wildcards: f"d{wildcards.d}"
    shell:
        "python src/mgatk_to_vireo.py {params.indir} {params.outdir} {params.donor}"


rule donor_copy_cells_meta:
    input: "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb"
    output:
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/donor{d}/mgatk_donor/cells_meta.tsv"
    params:
        cells_meta = lambda wildcards, input: join(dirname(input[0]), "cells_meta.tsv")
    shell: "cp {params.cells_meta} {output}"


rule vireo_mgatkdonor:
    input:
        mtx= "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/donor{d}/mgatk_donor/cellSNP.tag.AD.mtx",
        cells="{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/donor{d}/mgatk_donor/cells_meta.tsv",
    #output: "clones/variants_mgatkdonor/vireo/clones.ipynb"
    output:
        #"{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/clones.ipynb",
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}.labels.png",#, category="lineage", subcategory="donorMGATK"),
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}.variants.labels.png",# category="lineage", subcategory="donorMGATK"),
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}_cells_meta.tsv",
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        #N_DONORS =config["N_DONORS"], #config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "2b_MT_Lineage_Construct_mgatkDonors.ipynb"),
        d = lambda wildcards: wildcards.d,
        #workdir = os.getcwd(),
        nclones= lambda wildcards: wildcards.nclones
    threads: 8
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p donor {params.d} -p nclones {params.nclones} {params.notebook} {params.OUTDIR}/output.ipynb"


rule vireo_mgatkdonor_concat:
    input:
        expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/nclones{{nclones}}/donor{d}_cells_meta.tsv",
               d = np.arange(config["N_DONORS"])),
    output:
         "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/nclones{nclones}/cells_meta.tsv",
    run:
        all = []
        for i in input:
            all.append(pd.read_csv(i, sep='\t'))
        all = pd.concat(all, ignore_index=True).sort_values(["donor", "lineage", "donor_index", "lineage_index"])
        if 'level_0' in all.columns:
                all = all.drop('level_0', axis=1)
        ic('all')
        ic(all.head())

        all.to_csv(output[0], sep='\t', index=False)



rule multiplex_report:
    input:
        multiext("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/", "multiplex_AF_SNPs_all_afFilt.png", "multiplex_clusters_all.labels.png")
    output:
        report("{output}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/multiplex.pdf", category="Multiplex")
    shell:
        "convert {input} {output}"


rule donor_report:
    input:
        expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/dendrograms/figures/donor{{d}}_dendrogram.{suf}",
               suf=["png", "depth.png", "withHigh.png"]),
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/mgatk_donor/d{d}.variantQC.png",
        expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/donor{{d}}.labels.png",
               nclones=nclonelist, variants=params_clones["variants"]),
        expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/donor{{d}}.variants.labels.png",
               nclones=nclonelist, variants=params_clones["variants"])
    output:
        report("{output}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/donor{d}.pdf", category="Donor-specific")
    shell:
        "convert {input} {output}"


########################################################################
## Break up variants by clones and get the types of variants
rule clones_type_variants:
    input:
        "{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{{variants}}/vireo/nclones{{nclones}}/donor{d}.labels.png"
    output: "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/nclones{nclones}/variants.ipynb"
    params:
        notebook=join("src", "vireo", join("6_MT_Clones_variantTypes.ipynb")),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS = config["N_DONORS"],
        sample_names = ','.join(samples.index), # make it as a list
        nclones = lambda  wildcards: wildcards.nclones, #lambda wildcards: wildcards.nclones, #config['multiplex']["n_clone_list"],#lambda wildcards: wildcards.nclones,
        var_thresh=0.001,
        vars_to_plot=10
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR}  -p nclones {params.nclones} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} -p var_thresh {params.var_thresh} -p vars_to_plot {params.vars_to_plot} {params.notebook} {output} && jupyter nbconvert --to pdf {output}"


########################################################################
rule enrichment_vireo:
    version: '1.0' # no +1 norm error
    input:
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/cells_meta.tsv"
    output:
        report("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/enrichment/volcano_Fisher_foldNorm.png", category="enrichment"),
    params:
        clones_indir = lambda wildcards, input: dirname(input[0]),#lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        nclones = lambda wildcards: wildcards.nclones, #clones_cfg['vireo']["params"]["nclonelist"],
        script = join("src", "lineage", "lineage_enrichment.py"),
        samples=",".join(samples.index)
    shell: "python {params.script} {params.clones_indir} {params.OUTDIR} {params.nclones} {params.samples}"


# Bring enrichment results into same space
rule vireo_process:
    input:
        #expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{{variants}}/vireo/nclones{nclones}/{f_out}",
                    #zip,
                    #nclones=nclonelist, f_out=["enrichment/volcano_Fisher_foldNorm.png", "variants.ipynb"]),
        expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{{variants}}/vireo/nclones{nclones}/enrichment/volcano_Fisher_foldNorm.png",
                    nclones=nclonelist)
        #expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{{variants}}/vireo/nclones{nclones}/
    output:
        temp("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/temp/.tmp")

    shell: "touch {output}"


def get_cells_meta(wildcards):
    wc = wildcards
    return f"{wc.output}/data/merged/MT/cellr_{wc.cellrbc}/numread_{wc.num_read}/filters/minC{wc.mincells}_minR{wc.minreads}_topN{wc.topN}_hetT{wc.hetthresh}_hetC{wc.minhetcells}_hetCount{wc.hetcountthresh}_bq{wc.bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/donor{wc.d}/mgatk_donor/cells_meta.tsv"
    # if wildcards.variants == "mgatkdonor":
    #     return f"{wc.output}/data/merged/MT/cellr_{wc.cellrbc}/numread_{wc.num_read}/filters/minC{wc.mincells}_minR{wc.minreads}_topN{wc.topN}_hetT{wc.hetthresh}_hetC{wc.minhetcells}_hetCount{wc.hetcountthresh}_bq{wc.bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/donor{wc.d}/mgatk_donor/cells_meta.tsv"
    # elif wildcards.variants == "simple":
    #     return f"{wc.output}/data/merged/MT/cellr_{wc.cellrbc}/numread_{wc.num_read}/filters/minC{wc.mincells}_minR{wc.minreads}_topN{wc.topN}_hetT{wc.hetthresh}_hetC{wc.minhetcells}_hetCount{wc.hetcountthresh}_bq{wc.bqthresh}/mgatk/vireoIn/multiplex/cells_meta.tsv"

rule knn:
    input:
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/mgatk_donor/d{d}.variant.rds"
    output:
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/cell_clusters.tsv",
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/donor{d}.variants.labels.png"
    params:
        mgatk_in = lambda wildcards, input: input[0].replace(".variant.rds", ".af.tsv"),
        name = lambda wildcards: f"donor{wildcards.d}", #"/data2/mito_lineage/data/processed/mttrace/TcellDupi_may17_2021/MTblacklist/pre/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/pre.variant.rds"
        outdir = lambda wildcards, output: dirname(output[0]),
        kparam = lambda wildcards: wildcards.kparam,
        note = join(ROOT_DIR, "R_scripts", "knn_clones.ipynb"),
        #cells_meta = get_cells_meta,
        #d = lambda wildcards: wildcards.d,
    shell:
        "papermill -k ir -p mgatk_in {params.mgatk_in} -p name {params.name} -p outdir {params.outdir} -p kparam {params.kparam} {params.note} {params.outdir}/output.ipynb"


rule knn_process:
    input:
        expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/cell_clusters.tsv",
               kparam=params_clones["knn"]["params"]["resolution"], d=np.arange(config["N_DONORS"]))
    output: "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/temp/.tmp"#.pipeline"
    shell:
        "touch {output}"
    #outdir <- ""#"/data2/mito_lineage/Analysis/annotation/output/data/TcellDupi_may17_2021/MTblacklist/"
        # Cluster parameters
        #

# Sort of like an all rule to bring downstream results together.
rule complete_lineage:
    input:
        "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/temp/.tmp",
    output: "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/.completed"
    shell: "touch {output}"



rule terminate:
    input:
         rules.complete_lineage.output
    output:
          "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/.pipeline"
    shell: "touch {output}"