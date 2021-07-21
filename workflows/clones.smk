from os.path import join, dirname
import copy
import pandas as pd
import os
import numpy as np
import snakemake
#from snakemake import rules
wildcard_constraints:
    cellr='True|False'

# #####################################
# # Flatten config dictionary parameters
# if "mttrace" in config:
#     for c in config['mttrace']:
#         config[c] = copy.deepcopy(config['mttrace'][c])
# # if "filters" in config:
# #     for c in config['filters']:
# #         config[c] = copy.deepcopy(config['filters'][c])
#
# # Flatten dictionary parameters
# if "mgatk" in config:
#     for c in config['mgatk']:
#         config[c] = copy.deepcopy(config['mgatk'][c])
#####################################

#####################################
## Setup config variables
print('config')
print(config)
ref_fa = config["ref_fa"]
samples = pd.read_table(config["samples"], dtype=str,sep=',').set_index(["sample_name"], drop=False)
num_cells = config['multiplex']["pseudo_multiplex"]["num_cells"]
is_prop = config['multiplex']["pseudo_multiplex"]["is_proportional"]
res = config["results"]
ft = config["filters"]
#####################################


rule all:
    input:
        expand("{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx",
               sample=samples['sample_name'].values,
               results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh']),
        expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/enrichment/.status",
               results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh']),
        expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb",
                results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], n_clones=config['multiplex']["n_clone_list"]),
        expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/variants/variants.ipynb",
               results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh']),
        expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/dendrograms/af_dendro.ipynb",
               results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], n_clones=config['multiplex']["n_clone_list"]),
        expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/cells_BC.csv",
                results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], n_clones=config['multiplex']["n_clone_list"]),
        expand("{results}/{sample}/clones/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones.ipynb",
               sample=samples['sample_name'].values,
               results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh']),
        # expand("{results}/{sample}/clones/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/n_clones_{n_clones}/variants.ipynb",
        #        sample=samples['sample_name'].values, n_clones=config['multiplex']["n_clone_list"],
        #        results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
        #        het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'])
def get_input(wildcards, config, type='filters'):
    print(samples.loc[wildcards.sample, 'sample'])
    if type == 'mttrace':
        return config['files']['mttrace']['coverage'][samples.loc[wildcards.sample, 'sample']]
    if type == 'filters':
        print('here')
        print(config['files']['filters'])
        return config['files']['filters']['coverage'][samples.loc[wildcards.sample,'sample']]

from os.path import dirname
################################################################
## Import from prior snakefile which preprocesses the 10x output.
## Here, we redefine the input to be based on our config['files'] dictionary
from snakemake.utils import min_version
min_version("6.0")
module main_wf:
    snakefile: "./Snakefile"
    config: config


# use rule create_filters from main_wf as cf with:
#     input:
#         lambda wc: get_input(wc, config, 'mttrace')
#     output:
#         af_f = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/af_by_cell.tsv",
#         cov = "{results}/{sample}/filters//minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt"

use rule get_refAllele from main_wf with:
    output: "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/chrM_refAllele.txt"


def get_aggr(in_d):
    #bc = join(in_d, 'aggregate', 'outs', 'filtered_peak_bc_matrix', 'barcodes.tsv')
    aggr_csv = join(in_d, 'aggr.csv')
    return aggr_csv



rule create_filters:
    input:
        concat_dir = lambda wc: get_input(wc, config, 'mttrace')#"{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.coverage.strands.txt.gz"
    output:
        af_f = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/af_by_cell.tsv",
        cov = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt"
    params:
        concat_d = lambda wildcards, input: dirname(input.concat_dir),
        ref_fa = config['mt_ref_fa'],
        name = lambda wildcards: wildcards.sample,
        filt_params = lambda wc: f"{wc.min_cells} {wc.min_reads} {wc.topN} {wc.het_thresh} {wc.min_het_cells} {wc.het_count_thresh} {wc.bq_thresh}"
    resources:
        mem_mb=90000
    #log: "{results}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/.outcfg"
    shell: "python src/calculate_AF_by_cell.py {params.concat_d} {output.af_f} {params.ref_fa} {params.name} {params.filt_params}"# --log {log}"



use rule mgatk from main_wf as mgatk with:
    input:
        all = rules.create_filters.output.cov,# "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt",
        refAllele=rules.get_refAllele.output[0]
        #refAllele = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/chrM_refAllele.txt"
    output:
        "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.variant.rds"



# use rule mgatk_vireoIn from main_wf as mgatk_vireoIn with:
#     input:
#         "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.variant.rds"
#     output:
#         "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx"

# use rule merged from main_wf with:
#     input:
#         lambda wildcards: expand("{{results}}/{sample}/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx",
#                                         sample=samples["sample_name"].values) # , minC=wildcards.minC, minAF=wildcards.minAF)
#     output: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx"
#
# use rule multiplex from main_wf as multiplex with:
#     input:
#         "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx"
#     output:
#         "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/multiplex.ipynb"
################################################################


def get_filt(w):
    return w.min_cells, w.min_reads, w.topN, w.het_thresh, w.min_het_cells, w.het_count_thresh, w.bq_thresh



#
# rule get_refAllele:
#     #input: config["mt_ref_fa"],
#     params: config["chrM_refAllele"]
#     output: "{results}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/chrM_refAllele.txt"
#     shell: 'cp {params} {output}'
#
#
# rule mgatk:
#     """ Run both toSeurat and call variants in one script"""
#     input:
#         all = "{results}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt",
#         refAllele = "{results}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/chrM_refAllele.txt"
#     output: "{results}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.variant.rds"
#     params:
#         data_dir=lambda wildcards, input: dirname(input.all),
#         sample = lambda wildcards: wildcards.sample,
#     shell:
#         "./R_scripts/wrap_mgatk.R {params.data_dir} filter_mgatk/{params.sample} FALSE"
#
#
rule mgatk_to_vireoIn:
    input: "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.variant.rds"
    output: "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx"
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        sample = lambda wildcards: wildcards.sample
    shell:
        "python src/mgatk_to_vireo.py {params.indir} {params.outdir} {params.sample}"

#
rule merged:
    input:
        lambda wildcards: expand("{{results}}/{sample}/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx",
                                        sample=samples["sample_name"].values)# , minC=wildcards.minC, minAF=wildcards.minAF)
    #input: lambda wildcards: expand("data/{{prefix}}/chrM/{name}_cellSNP_minC{{mt_minC}}_minAF{{mt_minAF}}", name=config["samples"])
    params:
        num_cells=num_cells,
        is_prop = is_prop,
        indirs = lambda wildcards, input: [dirname(x) for x in input],
        outdir = lambda wildcards, output: dirname(output[0]),
        prefix= ','.join(samples["sample_name"])
    #log: "logs/{prefix}/vireo/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}.log"
    output: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx"
    shell: "python -m src.pseudo_batch {params.outdir} {params.indirs} --num_cells {params.num_cells} --is_prop {params.is_prop} --samples {params.prefix}" # > {log} 2>&1"


rule multiplex:
    input: rules.merged.output
    output:
        "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/multiplex.ipynb"
    params:
        N_DONORS=config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "1_MT_Donors_multiplex.ipynb" ),
        sample_names= ','.join(samples["sample_name"].values), # make it as a list
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        workdir = os.getcwd(),
        to_elbo = False
    shell: "papermill --cwd {params.workdir} -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p sample_names {params.sample_names} -p to_elbo {params.to_elbo} {params.notebook} {output}"


rule clones:
    input: rules.multiplex.output
    output: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/clones.ipynb"
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "2_MT_Lineage_Construct.ipynb"),
        workdir = os.getcwd()
    shell: "papermill --cwd {params.workdir} -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} {params.notebook} {output}"



rule enrichment:
    input:
        #rules.multiplex.output[0],
        rules.clones.output[0],
    output: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/enrichment/.status"
    params:
        clones_indir = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        n_clones = config['multiplex']["n_clone_list"],
        script = join("src", "vireo", "lineage_enrichment.py"),
        samples=",".join(config["multiplex"]['samples'])
    shell: "python {params.script} {params.clones_indir} {params.OUTDIR} {params.n_clones} {params.samples}"


rule donors_plotAF:
    input:
        #clone="data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/lineages/lineage.ipynb",
        mult=rules.multiplex.output[0]#"data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/multiplex.ipynb"
    output:
        "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/dendrograms/af_dendro.ipynb",
        #report("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/dendrograms/figures/")
        report(expand("{{results}}/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/dendrograms/figures/donor{n}_dendrogram.png",
               n=np.arange(config["multiplex"]["N_DONORS"])))
    params:
        #INDIR = lambda wildcards, input: dirname(input[0]),
        INDIR = lambda wildcards, input: dirname(input.mult),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config['multiplex']['N_DONORS'], #config["multiplex"]["N_DONORS"],
        sample_names= ",".join(config["multiplex"]['samples']), # make it as a list
        notebook=join("src", "vireo", "3_MT_Donors_Dendrogram.ipynb"),
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} {params.notebook} {output[0]}"


rule donors_type_variants:
    input: rules.multiplex.output[0]
    params:
        notebook=join("src", "vireo", join("5_MT_Donors_variantTypes.ipynb")),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["multiplex"]["N_DONORS"],
        sample_names= ','.join(samples["sample_name"].values) # make it as a list
    output: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/variants/variants.ipynb",
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} {params.notebook} {output} && jupyter nbconvert --to pdf {output}"


rule clones_plotAF:
    input:
        #clone="data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/lineages/lineage.ipynb",
        mult = rules.multiplex.output[0]#"data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/multiplex.ipynb"
    output:
        "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/dendrograms/dendrograms.ipynb",
        #report("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/nfigures/donor{n}_dendrogram.png")
        report(expand("{{results}}/data/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/dendrograms/figures/donor{n}_dendrogram.png",
                n=np.arange(config["multiplex"]["N_DONORS"]),
                n_clones=config['multiplex']["n_clone_list"]))
    params:
        #INDIR = lambda wildcards, input: dirname(input[0]),
        INDIR = lambda wildcards, input: dirname(input.mult),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["multiplex"]["N_DONORS"],
        sample_names= ','.join(samples["sample_name"].values), # make it as a list
        notebook=join("src", "vireo", "3_MT_Donors_Dendrogram.ipynb"),
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} {params.notebook} {output[0]}"


rule clones_type_variants:
    input: rules.clones.output[0]
    params:
        notebook=join("src", "vireo", join("6_MT_Clones_variantTypes.ipynb")),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS = config["multiplex"]["N_DONORS"],
        sample_names = ','.join(config['multiplex']["samples"]), # make it as a list
        n_clones = lambda wildcards: wildcards.n_clones, #config['multiplex']["n_clone_list"],#lambda wildcards: wildcards.n_clones,
        var_thresh=0.001,
        vars_to_plot=10
    #output: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/variants/variants.ipynb",
    output: "{results}/{sample}/clones/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb"
    #output: report(lambda wildcards: expand("{{results}}/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb", n_clones=config['multiplex']["n_clone_list"]))
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR}  -p n_clones {params.n_clones} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} -p var_thresh {params.var_thresh} -p vars_to_plot {params.vars_to_plot} {params.notebook} {output} && jupyter nbconvert --to pdf {output}"



rule clones_exp:
    input: rules.multiplex.output
    output: "{results}/{sample}/clones/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones.ipynb"
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "2_MT_Lineage_Construct_separateCond.ipynb"),
        workdir = os.getcwd(),
        exp = lambda wildcards: wildcards.sample
    shell: "papermill --cwd {params.workdir} -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p exp {params.exp} -p N_DONORS {params.N_DONORS} {params.notebook} {output} && jupyter nbconvert --to pdf {output}"


rule clones_type_variants_exp:
    input: rules.clones_exp.output[0]
    params:
        notebook=join("src", "vireo", join("6_MT_Clones_variantTypes_separateCond.ipynb")),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS = config["multiplex"]["N_DONORS"],
        sample_names = ','.join(config['multiplex']["samples"]), # make it as a list
        n_clones = lambda wildcards: wildcards.n_clones, #config['multiplex']["n_clone_list"],#lambda wildcards: wildcards.n_clones,
        var_thresh=0.001,
        vars_to_plot=10,
        exp = lambda wildcards: wildcards.sample
    #output: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/variants/variants.ipynb",
    output: "{results}/{sample}/clones/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/n_clones_{n_clones}/variants.ipynb"
    #output: report(lambda wildcards: expand("{{results}}/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb", n_clones=config['multiplex']["n_clone_list"]))
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p n_clones {params.n_clones} -p exp {params.exp} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} -p var_thresh {params.var_thresh} -p vars_to_plot {params.vars_to_plot} {params.notebook} {output} && jupyter nbconvert --to pdf {output}"


use rule merge_lineage_nuclear_barcodes from main_wf with:
    input:
        cl="{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/clones.ipynb",
        aggr=get_aggr(config['mtscATAC_OUTDIR']) #lambda wildcards: join(config['mtscATAC_OUTDIR'], 'aggr.csv')
    output:
        bc="{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/cells_BC.csv",
        meta="{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/cells_meta_aggr.tsv"

# rule merge_lineage_nuclear_barcodes:
#     input:
#         "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/clones.ipynb",
#         get_aggr(config['mtscATAC_OUTDIR'])
#     params:
#         N_DONORS=config["multiplex"]["N_DONORS"],
#         cells_meta = lambda wildcards, output: join(dirname(input[0]),f'lineage{wildcards.n_clones}', "cells_meta.tsv"),
#     output:
#         "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/cells_BC.csv"
#         "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/cells_meta_aggr.tsv"
#     run:
#         cells_meta = pd.read_csv(params.cells_meta)
#         aggr = pd.read_csv(input[1])
#         aggr.index = aggr.index + 1 #1-based index
#         samples_d = {val:i for i, val in enumerate(aggr['library'].values)}
#         cells_meta['BC aggr'] = cells_meta.apply(lambda x: f"{x['raw ID'].replace('-1','')}-{samples_d[x['condition']]}", axis=1)
#         cells_meta.to_csv(output[1], sep='\t', index=False)
#         cells_meta['BC aggr'].to_csv(output[0], index=False, header=False)
