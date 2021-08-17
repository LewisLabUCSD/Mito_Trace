import os
import pandas as pd
from snakemake.utils import min_version
min_version("6.0")
from os.path import dirname, join
import numpy as np
# module main_wf:
#     snakefile: "./Snakefile"
#     config: config

## Setup config variables
print('config')
print(config)
ref_fa = config["ref_fa"]
samples = config['multiplex']["samples"] #pd.read_table(config["samples"], dtype=str,sep=',').set_index(["sample_name"], drop=False)
num_cells = config['multiplex']["pseudo_multiplex"]["num_cells"]
is_prop = config['multiplex']["pseudo_multiplex"]["is_proportional"]
res = config["results"]

rule all:
    input:
        expand("{results}/merged/enrichment/.status", results=res),
        expand("{results}/merged/clones/n_clones_{n_clones}/variants.ipynb",
                results=res),
        expand("{results}/merged/multiplex/variants/variants.ipynb",
               results=res),
        expand("{results}/merged/dendrograms/af_dendro.ipynb",
               results=res, n_clones=config['multiplex']["n_clone_list"]),
        expand("{results}/merged/clones/n_clones_{n_clones}/cells_BC.csv",
                results=res),
        expand("{results}/{sample}/clones/clones.ipynb",
               sample=samples['sample_name'].values,
               results=res),

rule merged:
    input:
        lambda wildcards: expand("{{results}}/{sample}/cellSNP.tag.AD.mtx",
                                 sample=samples["sample_name"].values)# , minC=wildcards.minC, minAF=wildcards.minAF)
    #input: lambda wildcards: expand("data/{{prefix}}/chrM/{name}_cellSNP_minC{{mt_minC}}_minAF{{mt_minAF}}", name=config["samples"])
    params:
        num_cells=num_cells,
        is_prop = is_prop,
        indirs = lambda wildcards, input: [dirname(x) for x in input],
        outdir = lambda wildcards, output: dirname(output[0]),
        prefix= ','.join(samples["sample_name"])
    #log: "logs/{prefix}/vireo/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}.log"
    output: "{results}/merged/cellSNP.tag.AD.mtx"
    shell: "python -m src.pseudo_batch {params.outdir} {params.indirs} --num_cells {params.num_cells} --is_prop {params.is_prop} --samples {params.prefix}" # > {log} 2>&1"


rule multiplex:
    input: rules.merged.output
    output:
        "{results}/merged/multiplex/multiplex.ipynb"
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
    output: "{results}/merged/clones/clones.ipynb"
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
    output: "{results}/merged/enrichment/.status"
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
        "{results}/merged/dendrograms/af_dendro.ipynb",
        #report("{results}/merged/dendrograms/figures/")
        report(expand("{{results}}/merged/dendrograms/figures/donor{n}_dendrogram.png",
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
    output: "{results}/merged/multiplex/variants/variants.ipynb",
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} {params.notebook} {output} && jupyter nbconvert --to pdf {output}"


rule clones_plotAF:
    input:
        #clone="data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/lineages/lineage.ipynb",
        mult = rules.multiplex.output[0]#"data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/multiplex.ipynb"
    output:
        "{results}/merged/clones/dendrograms/dendrograms.ipynb",
        #report("{results}/merged/clones/nfigures/donor{n}_dendrogram.png")
        report(expand("{{results}}/data/merged/clones/n_clones_{n_clones}/dendrograms/figures/donor{n}_dendrogram.png",
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
    #output: "{results}/merged/clones/variants/variants.ipynb",
    output: "{results}/merged/clones/n_clones_{n_clones}/variants.ipynb"
    #output: report(lambda wildcards: expand("{{results}}/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb", n_clones=config['multiplex']["n_clone_list"]))
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR}  -p n_clones {params.n_clones} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} -p var_thresh {params.var_thresh} -p vars_to_plot {params.vars_to_plot} {params.notebook} {output} && jupyter nbconvert --to pdf {output}"
