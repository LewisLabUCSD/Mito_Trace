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
if "filters" in config:
    ft = config["filters"]
else:
    ft = config["mttrace"]
#####################################


rule all:
    input:
        expand("{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx",
               sample=samples['sample_name'].values,
               results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh']),
        # expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/enrichment/.status",
        #        results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
        #        het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh']),
        expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/enrichment/volcano_donor{d}.clones{n_clone}_Fisher_foldNorm.png",
                d=range(config["multiplex"]["N_DONORS"]), n_clone=config['multiplex']['n_clone_list'],
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
         expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/donor_mgatk/vireoIn/clones/clones.ipynb",
                results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
                het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], d=range(config["multiplex"]["N_DONORS"])),
         expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/dendrograms/figures/donor{n}_dendrogram.png",
                n=np.arange(config["multiplex"]["N_DONORS"]), results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
                het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], d=range(config["multiplex"]["N_DONORS"]))

        # expand("{results}/clones_knn/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/clones_knn.ipynb",
        #        results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
        #        het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], d=np.arange(config["multiplex"]["N_DONORS"]))
        # expand("{results}/{sample}/clones/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/n_clones_{n_clones}/variants.ipynb",
        #        sample=samples['sample_name'].values, n_clones=config['multiplex']["n_clone_list"],
        #        results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
        #        het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'])

def get_input(wildcards, config, type='filters'):
    print(samples.loc[wildcards.sample, 'sample'])
    if type == 'mttrace':
        if samples.loc[wildcards.sample, 'sample_name'] in  config['files']['mttrace']['coverage']:
            return config['files']['mttrace']['coverage'][samples.loc[wildcards.sample, 'sample_name']]
        else:
            return config['files']['mttrace']['coverage'][samples.loc[wildcards.sample, 'sample']]
    if type == 'filters':
        print('here')
        print(config['files']['filters'])
        if samples.loc[wildcards.sample, 'sample_name'] in  config['files']['filters']['coverage']:
            return config['files']['filters']['coverage'][samples.loc[wildcards.sample,'sample_name']]
        else:
            return config['files']['filters']['coverage'][samples.loc[wildcards.sample, 'sample']]

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
        #af_f = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/af_by_cell.tsv",
        cov = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt",
        stats = report(multiext("{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/", "stats.csv", "initial_cell_depth.png","heatmap.png"))
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
        vars = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.variant.rds",
        vars_qc = report("{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.variantQC.png")



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
        "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/multiplex.ipynb",
        results = report(multiext("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/", "multiplex_AF_SNPs_all_afFilt.png", "multiplex_clusters_all.labels.png"))
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
    input: rules.multiplex.output[0]
    output:
        #"{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/clones.ipynb",
        report(expand("{{results}}/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/clones/lineage{n_clone}/donor{d}_lineage_{n_clone}OUT.labels.png",
                        d=range(config["multiplex"]["N_DONORS"]), n_clone=config['multiplex']['n_clone_list'])),
        report(expand("{{results}}/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/clones/lineage{n_clone}/donor{d}_lineage_{n_clone}OUT.variants.labels.png",
               d=range(config["multiplex"]["N_DONORS"]), n_clone=config['multiplex']['n_clone_list']))
    params:
        output_notebook = lambda wildcards, output: join(dirname(dirname(output[0])), 'clones.ipynb'),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "2_MT_Lineage_Construct.ipynb"),
        workdir = os.getcwd(),
        n_clone= ",".join([str(x) for x in config['multiplex']['n_clone_list']])
    shell: "papermill --cwd {params.workdir} -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p n_clone_list {params.n_clone} {params.notebook} {params.output_notebook}"



rule enrichment:
    input:
        #rules.multiplex.output[0],
        #rules.clones.output[0],
        expand("{{results}}/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/clones/lineage{n_clone}/donor{d}_lineage_{n_clone}OUT.labels.png",
            d=range(config["multiplex"]["N_DONORS"]), n_clone=config['multiplex']['n_clone_list'])
    output:
        #"{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/enrichment/.status",
        report(expand("{{results}}/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/enrichment/volcano.clones{n_clone}_Fisher_foldNorm.png",
                        n_clone=config['multiplex']['n_clone_list'])),
        # report(expand("{{results}}/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/enrichment/enrichment_donor{n}_clones{n_clone}_FisherNorm.csv",
        #                 d=range(config["multiplex"]["N_DONORS"]), n_clone=",".join([str(x) for x in config['multiplex']['n_clone_list']])))

        #report("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/enrichment/enrichment_volcano_donor{n}.clones{n_clone}_FisherNorm.png"),
        #report("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/enrichment/enrichment_donor{n}_clones{n_clone}_FisherNorm.csv")
    params:
        clones_indir = lambda wildcards, input: dirname(dirname(dirname(input[0]))),#lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        n_clones = config['multiplex']["n_clone_list"],
        script = join("src", "lineage_enrichment.py"),
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
        report(expand("{{results}}/data/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/dendrograms/figures/donor{d}_dendrogram.png",
                d=np.arange(config["multiplex"]["N_DONORS"]),
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



########################################################################
## Workflow B: After multiplexing, separate by donor, grab variants from filters
## that overlap with both conditions, and then run mgatk to call variants again.
##
# Extract pileups for each donor from original filter matrices
rule scPileup_filter_mgatk:
    input:
        cov = expand("{{results}}/{sample}/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/{sample}.coverage.txt", sample=samples['sample_name'].values),
        #cov = expand("{{results}}/data/{sample}/MT/cellr_{{cellr_bc}}/{sample}_{{num_read}}/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/{sample}.coverage.txt", sample=samples['sample_name'].values),
        mult = "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/multiplex.ipynb"
    output:
        cov = expand("{{results}}/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/d{d}.coverage.txt",
                     d=range(config["multiplex"]["N_DONORS"]))
    params:
        #outdir = lambda wildcards, output: dirname(output.cov),
        cells_meta = lambda wildcards, input: join(dirname(input.mult), "cells_meta.tsv"),
        sample = samples['sample_name'].values

    run:
        cells_meta = pd.read_csv(params.cells_meta, sep='\t')
        # Save pileup for each donor. Loop through samples
        # and get overlapping positions
        for d, df in cells_meta.groupby('donor'):
            d = int(d)
            nt_pileups = {}
            cond_positions = {}
            nt_pileups["coverage"] = pd.DataFrame()

            for ind, curr_cov_f in enumerate(input.cov):
                print('curr_cov_f ', curr_cov_f )
                curr_samp_cov = pd.read_csv(curr_cov_f, header=None) #load sample coverage
                condition = params.sample[ind]
                curr_samp_donor_df = df[df["condition"]==condition] #Condition-donors
                filt_df = curr_samp_cov.loc[curr_samp_cov[1].isin(curr_samp_donor_df["raw ID"])]
                #print('filt_df', filt_df.head())
                filt_df[1] = filt_df[1] + "_" + params.sample[ind]
                cond_positions[params.sample[ind]] = set(filt_df[0].values)
                nt_pileups["coverage"] = nt_pileups["coverage"].append(filt_df, ignore_index=True)

            # Filter by overlapping positions:
            pos_to_keep = list(cond_positions.values())[0].intersection(*list(cond_positions.values()))
            nt_pileups["coverage"] = nt_pileups["coverage"][nt_pileups["coverage"][0].isin(pos_to_keep)]

            # Some cells may not be present anymore, need to store them.
            cells_to_keep = set(nt_pileups["coverage"][1])

            # Save coverage and depth
            nt_pileups["coverage"].to_csv(str(output[int(d)]), index=False,
                                          header=False)
            depth = nt_pileups["coverage"].groupby(1).sum()[2]/16569
            depth.to_csv(join(dirname(output[int(d)]), f"d{d}.depthTable.txt"), sep='\t', header=False)

            # Filter for the nucleotides as well
            for nuc in ["C", "A", "G", "T"]:
                nt_pileups[nuc] = pd.DataFrame()
                for ind, curr_cov_f in enumerate(input.cov):
                    curr_samp_cov = pd.read_csv(join(dirname(curr_cov_f), f"{params.sample[ind]}.{nuc}.txt"), header=None) #load sample coverage
                    condition = params.sample[ind]
                    curr_samp_donor_df = df[df["condition"]==condition]
                    filt_df = curr_samp_cov.loc[curr_samp_cov[1].isin(curr_samp_donor_df["raw ID"])] # Filter by donor
                    ## Filter by positions to keep and cells
                    filt_df[1] = filt_df[1] + "_" + params.sample[ind]
                    filt_df = filt_df[((filt_df[0].isin(pos_to_keep)) & (filt_df[1].isin(cells_to_keep)))]
                    nt_pileups[nuc] = nt_pileups[nuc].append(filt_df, ignore_index=True)
                # Save coverage
                nt_pileups[nuc].to_csv(join(dirname(output[int(d)]), f"d{d}.{nuc}.txt"), header=False, index=False)

            cells_meta = cells_meta[cells_meta["ID"].isin(cells_to_keep)]
            cells_meta.to_csv(join(dirname(output[int(d)]), "cells_meta.tsv"), sep='\t', index=False)
    #shell: "python donor_filter_mgatk.py {params.indir} {params.cells_meta} {params.outdir}"


rule donor_get_refAllele:
    #input: config["mt_ref_fa"],
    params: config["chrM_refAllele"]
    output: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/chrM_refAllele.txt"
    shell: 'cp {params} {output}'


rule donor_mgatk:
    input:
        "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/d{d}.coverage.txt",
        rules.donor_get_refAllele.output
    output: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/donor_mgatk/d{d}.variant.rds"
    params:
        data_dir =lambda wildcards, input: dirname(input[0]),
        donor = lambda wildcards: f"d{wildcards.d}",
    shell:
        "./R_scripts/wrap_mgatk.R {params.data_dir} donor_mgatk/{params.donor} FALSE"


rule donor_mgatk_to_vireoIn:
    input: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/donor_mgatk/d{d}.variant.rds"
    output: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/donor_mgatk/vireoIn/cellSNP.tag.AD.mtx"
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        donor = lambda wildcards: f"d{wildcards.d}"
    shell:
        "python src/mgatk_to_vireo.py {params.indir} {params.outdir} {params.donor}"

rule donor_copy_cells_meta:
    input: rules.multiplex.output
    output:
        "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/donor_mgatk/vireoIn/cells_meta.tsv"
    params:
        cells_meta = lambda wildcards, input: join(dirname(input[0]), "cells_meta.tsv")
    shell: "cp {params.cells_meta} {output}"


rule clones_donor_mgatk:
    input:
        "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/donor_mgatk/vireoIn/cellSNP.tag.AD.mtx",
        "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/donor_mgatk/vireoIn/cells_meta.tsv"
    output: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/donor_mgatk/vireoIn/clones/clones.ipynb"
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        donor=lambda wildcards: wildcards.d, #config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "2b_MT_Lineage_Construct_mgatkDonors.ipynb"),
        workdir = os.getcwd(),
        n_clone=",".join([str(x) for x in config['multiplex']['n_clone_list']])
    shell: "papermill --cwd {params.workdir} -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p donor {params.donor} -p n_clone_list {params.n_clone} {params.notebook} {output}"


## Workflow B ends here.
########################################################################
########################################################################


########################################################################
## Workflow C: Run lineage tracing on each sample separately after multiplex
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
########################################################################
########################################################################

########################################################################
## Workflow D: Use knn clustering in seurat to detect clones
rule merge_mgatk:
    input:
         #"{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.variant.rds"
        expand("{{results}}/{sample}/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/{sample}.variant.rds",
               sample=samples["sample_name"].values)
    output:
        mgatk="{results}/clones_knn/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk.variants.merged.rds"
    shell: "Rscript ./R_scripts/mergeVariants.R {input} {output.mgatk}"


rule clones_knn:
    input:
        "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/multiplex.ipynb",
        mgatk="{results}/clones_knn/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk.variants.merged.rds"
    output:
        "{results}/clones_knn/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/clones_knn.ipynb"
    params:
        cells_f = lambda wildcards, input: join(dirname(input.mgatk),"cells_meta.tsv"),
        outdir = lambda wildcards, output: dirname(output[0]),
        notebook = join("R_scripts", "call_clones.ipynb"),
        cells_col = "donor"
    shell:
        "papermill -p mgatk_in {input.mgatk} -p outdir {params.outdir} -p cells_f {params.cells_f} -p cells_col {params.cells_col} {params.notebook} {output}  && jupyter nbconvert --to pdf {output}"
########################################################################
########################################################################



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
