wildcard_constraints:
    mapq="\d+",
    cellr='True|False'

from os.path import join, dirname
import pandas as pd
import copy
from src.utils.paramspace import Paramspace
from src.config import ROOT_DIR

outputs = {}
for ind,p in enumerate((config["pipeline"])):
    f"{ind:.2d}_{config[p]}"
    outputs[p] = Paramspace(config[p], param_sep="_" )

pipeline = config["pipeline"]
#for p in pipeline:


########
rule all:
    input:
        expand(outputs)
#         expand("{results}/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz",
#                results=res,sample=samples["sample_name"].values, num_read=config["num_reads_filter"]),
#         expand("{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/mgatk/{sample}.lowC{low_cov_thresh}.variant.rds",
#                results=res,sample=samples["sample_name"].values,
#                cellr_bc=cellr_bc, num_read=num_reads_filter, low_cov_thresh=config["low_cov_thresh"]),
# #         expand("{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}_MT_position_coverage.png",
# #                results=res,sample=samples["sample_name"].values, num_read=config["num_reads_filter"], cellr_bc=config["use_cellr_barcode"])
# # ########


rule filters:
    input:
        concat_dir = "{indir}/scPileup/{sample}.coverage.strands.txt.gz"
    #output:  "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/af_by_cell.tsv"
    output:
        af_f = "{indir}/filters/af_by_cell.tsv",
        cov = "{indir}/filters/{sample}.coverage.txt"
    params:
        concat_d = lambda wildcards, input: dirname(input.concat_dir),
        ref_fa = config['mt_ref_fa'],
        name = lambda wildcards: wildcards.sample,
        filt_params = get_filt,
    resources:
        mem_mb=90000
    log: "{indir}/logs/filters/.log"
    shell: "python src/calculate_AF_by_cell.py {params.concat_d} {output.af_f} {params.ref_fa} {params.name} {params.filt_params} --log {log}"

import numpy as np
def params_string(name, params, name_sep="__",params_sep="_"):
    curr = f"{name}{name_sep}"
    for p in params:
        curr = f"{curr}{p}{params[p]}{params_sep}"
    curr = curr.strip(params_sep)
    return curr

def get_input(config, type):
    pipe = config['pipe']
    ind = np.flatnonzero(pipe==type)
    if ind == 0:
        return rules.filters.output[0]
    else:
        return f"{params_string(pipe[ind-1], config[pipe[ind-1]])}/.outcfg"
        #return f"{pipe[ind-1]}/.outcfg"


rule variants:
    """ Run both toSeurat and call variants in one script"""
    input:
        get_input(config['pipeline'], "variants")
        #all = "{indir}/filters/{sample}.coverage.txt", #get_input(variants)
        #refAllele = "{indir}/filters/chrM_refAllele.txt"
    output: "{indir}/mgatk/{sample}.variant.rds"
    params:
        data_dir=lambda wildcards, input: dirname(input.all),
        sample = lambda wildcards: wildcards.sample,
    shell:
        "./R_scripts/wrap_mgatk.R {params.data_dir} filter_mgatk/{params.sample} FALSE"

rule mgatk_to_vireoIn:
    input: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.variant.rds"
    output: directory("{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/")
    params:
        AF= lambda wildcards, input: join(dirname(input[0]), f'{wildcards.sample}.af.tsv')
    shell:
        "python src/mgatk_to_vireo.py {params.AF} {output}"



## ChrMT Pseudo-multiplex
rule pseudo:
    """Subsets and combines the samples to see if it can properly
        deconvolve the two."""
    #TODO
    input:
        #lambda wildcards: expand("data/{{prefix}}/chrM/{name}_cellSNP_minC{{mt_minC}}_minAF{{mt_minAF}}", name=config["samples"])# , minC=wildcards.minC, minAF=wildcards.minAF)
        lambda wildcards: expand("{indir}/{{prefix}}/chrM/{name}_cellSNP_minC{{mt_minC}}_minAF{{mt_minAF}}") #get_input(config['pipeline'], "pseudo")
    params:
        num_cells=config['num_cells'],
        is_prop = config['is_prop'],
        outdir = lambda wildcards, output: dirname(output[0]),
        prefix= ','.join(config["samples"]) # make it as a list
    output: "{indir}/multiplex/cellSNP.tag.AD.mtx"
    log: "logs/{prefix}/vireo/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}.log"
    shell: "python -m src.pseudo_batch {params.outdir} {input} --num_cells {params.num_cells} --is_prop {params.is_prop} --samples {params.prefix} > {log} 2>&1"


rule multiplex:
    input:
        get_input(config['pipeline'], "multiplex")
        #AD_F="{indir}/multiplex/cellSNP.tag.AD.mtx",
    output:
        out_note="{indir}/multiplex/multiplex.ipynb",
    params:
        N_DONORS=config["N_DONORS"],
        notebook=join(ROOT_DIR, "src", "vireo", "1_MT_Donors_multiplex.ipynb" ),
        sample_names= ','.join(config["samples"]), # make it as a list
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output.out_note),
        #sample_csv = config["sample_csv"]
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p sample_names {params.sample_names} {params.notebook} {output}"


rule clones:
    input: f"{config['pipeline'][-1]}/.outcfg"#"{indir}/multiplex/multiplex.ipynb"
    output: "{indir}/lineages/lineage.ipynb"
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["N_DONORS"],
        notebook=join(ROOT_DIR, "src", "vireo", "2_MT_Lineage_Construct.ipynb"),
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} {params.notebook} {output}"


rule summary:
    input:
        "{indir}/lineages/lineage.ipynb" #f"{config['pipeline'][-1]}/.outcfg"#mult="{indir}/multiplex/multiplex.ipynb"
    output: "{indir}/multiplex/dendrograms/af_dendro.ipynb"
    params:
        #INDIR = lambda wildcards, input: dirname(input[0]),
        INDIR = lambda wildcards, input: dirname(input.mult),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["N_DONORS"],
        sample_names= ','.join(config["samples"]), # make it as a list
        notebook=join(ROOT_DIR, "src", "vireo", "3_MT_Donors_Dendrogram.ipynb"),
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS}  {params.notebook} {output}"


