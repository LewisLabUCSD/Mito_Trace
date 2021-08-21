#from src.config import ROOT_DIR
#workdir: ROOT_DIR
wildcard_constraints:
    cellr='True|False'

from src.config import ROOT_DIR
from src.utils.parse_config import read_config_file
import os
from os.path import join, dirname
import pandas as pd
import copy
from snakemake.utils import min_version
min_version("6.0")

res = config["results"]
params = read_config_file(config["params"])
print('config')
print(config)

params_mt = params["mtpreproc"]
params_ft = params["filters"]
params_mult = params["multiplex"]

cellr_bc = params_mt["use_cellr_barcode"]
num_reads_filter = params_mt["num_reads_filter"]
maxBP = config["maxBP"]
ref_fa = config["ref_fa"]
samples = pd.read_table(config["samples"], dtype=str,sep=',').set_index(["sample_name"], drop=False)
ncells_thresh_mgatk=config['ncells_thresh_mgatk']
num_cells = config['multiplex']["pseudo_multiplex"]["num_cells"]
is_prop = config['multiplex']["pseudo_multiplex"]["is_proportional"]


#workdir: config["work_dir"]
rule all:
    input:
         # A. Preproc: Initial MT coverage per position per cell
        expand("{results}/data/{sample}/MT/cellr_{cellr_bc}/numread_{num_read}_MT_position_coverage.png",
               results=res,sample=samples["sample_name"].values, num_read=config["num_reads_filter"], cellr_bc=config["use_cellr_barcode"]),
         # B. Multiplex: Multiplexed VAF and depth dendrograms for each donor and selected variants
        expand("{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk/vireoIn/dendrograms/af_dendro.ipynb",
               results=res,cellr_bc=cellr_bc, num_read=num_reads_filter,
               min_cells=params_ft['min_cells'],min_reads=params_ft['min_reads'],topN=params_ft["topN"],het_thresh=params_ft['het_thresh'],min_het_cells=params_ft['min_het_cells'],
               het_count_thresh=params_ft['het_count_thresh'], bq_thresh=params_ft['bq_thresh']),
        # C. Clones: Variant types
        expand("{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb",
               results=res, cellr_bc=cellr_bc, num_read=num_reads_filter,
               min_cells=params_ft['min_cells'],min_reads=params_ft['min_reads'],topN=params_ft["topN"],het_thresh=params_ft['het_thresh'],min_het_cells=params_ft['min_het_cells'],
               het_count_thresh=params_ft['het_count_thresh'], bq_thresh=params_ft['bq_thresh'], n_clones=config['multiplex']['n_clone_list']),
        # D. Enrichment: Volcano plot
        expand("{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk/vireoIn/enrichment/.status",
               results=res, cellr_bc=cellr_bc, num_read=num_reads_filter,
               min_cells=params_ft['min_cells'],min_reads=params_ft['min_reads'],topN=params_ft["topN"],het_thresh=params_ft['het_thresh'],min_het_cells=params_ft['min_het_cells'],
               het_count_thresh=params_ft['het_count_thresh'], bq_thresh=params_ft['bq_thresh']),
        # E. Lineage labels for cells, which can be input to loupe browser
        expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/cells_BC.csv",
                results=res, min_cells=params_ft['min_cells'],min_reads=params_ft['min_reads'],topN=params_ft["topN"],het_thresh=params_ft['het_thresh'],min_het_cells=params_ft['min_het_cells'],
                het_count_thresh=params_ft['het_count_thresh'], bq_thresh=params_ft['bq_thresh'], n_clones=config['multiplex']["n_clone_list"]),

        ## Workflow B: Re-running mgatk after multiplexing with the intersection of the filtered variants (not the mgatk variants subset). This is to see if any mutation
                        # within a donor has a high variability, which might've been masked before. It also puts back the depth which was 0'd out in variants not seen in a sample, which could have biased the results.
         # F. Clones~mgatkRerun: types of clones
        expand("{results}}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/lineage{n_clone}/donor{d}_lineage_{n_clone}OUT.variants.labels.png",
               d=range(config["N_DONORS"]), n_clone=config['multiplex']['n_clone_list'],
               results=res, min_cells=params_ft['min_cells'],min_reads=params_ft['min_reads'],topN=params_ft["topN"],het_thresh=params_ft['het_thresh'],min_het_cells=params_ft['min_het_cells'],
               het_count_thresh=params_ft['het_count_thresh'], bq_thresh=params_ft['bq_thresh']
               ),
         # G. Clones~mgatkRerun: types of clones
        expand("{results}}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/lineage{n_clone}/donor{d}_lineage_{n_clone}OUT.variants.labels.png",
               d=range(config["N_DONORS"]), n_clone=config['multiplex']['n_clone_list'],
               results=res, min_cells=params_ft['min_cells'],min_reads=params_ft['min_reads'],topN=params_ft["topN"],het_thresh=params_ft['het_thresh'],min_het_cells=params_ft['min_het_cells'],
               het_count_thresh=params_ft['het_count_thresh'], bq_thresh=params_ft['bq_thresh']
               )

        # H. Clones~mgatkRerun: cl


################################################################
## Import from prior snakefile modules
## Here, we redefine the input to be based on our config['files'] dictionary
from snakemake.utils import min_version
min_version("6.0")
module mtpreprocMod:
    snakefile: "./mt_preprocess.Snakefile"
    config: config

module mgatkMod:
    snakefile: "./mgatk.Snakefile"
    config: config

module multMode:
    snakefile: "./multiplex.Snakefile"
    config: config

module lineageMod:
    snakefile: "./lineage.Snakefile"
    config: config


use rule all from mtpreprocMod as mgatk with:
    input:
        all = rules.create_filters.output.cov,# "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt",
        refAllele=rules.get_refAllele.output[0]
    output:
        vars = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk/{sample}.variant.rds",
        vars_qc = report("{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk/{sample}.variantQC.png")
    sample = lambda wildcards: wildcards.sample,


use rule all from mgatkMod with:
    input:
        all = rules.create_filters.output.cov,# "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt",
        refAllele=rules.get_refAllele.output[0]
        #refAllele = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/chrM_refAllele.txt"
    output:
        vars = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk/{sample}.variant.rds",
        vars_qc = report("{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk/{sample}.variantQC.png")


rule merged:
    """ Merge the sample pileup matrices
    """
    input:
        lambda wildcards: expand("{{results}}/{sample}/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/mgatk/vireoIn/cellSNP.tag.AD.mtx",
                                        sample=samples["sample_name"].values)# , minC=wildcards.minC, minAF=wildcards.minAF)
    #input: lambda wildcards: expand("data/{{prefix}}/chrM/{name}_cellSNP_minC{{mt_minC}}_minAF{{mt_minAF}}", name=config["samples"])
    params:
        num_cells=num_cells,
        is_prop = is_prop,
        indirs = lambda wildcards, input: [dirname(x) for x in input],
        outdir = lambda wildcards, output: dirname(output[0]),
        prefix= ','.join(samples["sample_name"])
    #log: "logs/{prefix}/vireo/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}.log"
    output: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk/vireoIn/cellSNP.tag.AD.mtx"
    shell: "python -m src.pseudo_batch {params.outdir} {params.indirs} --num_cells {params.num_cells} --is_prop {params.is_prop} --samples {params.prefix}" # > {log} 2>&1"



use rule all from multMod with:
    input:
        all = rules.create_filters.output.cov,# "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt",
        refAllele=rules.get_refAllele.output[0]
        #refAllele = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/chrM_refAllele.txt"
    output:
        vars = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk/{sample}.variant.rds",
        vars_qc = report("{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk/{sample}.variantQC.png")



use rule all from lineageMod with:
    input:
        all = rules.create_filters.output.cov,# "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt",
        refAllele=rules.get_refAllele.output[0]
        #refAllele = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/chrM_refAllele.txt"
    output:
        vars = "{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk/{sample}.variant.rds",
        vars_qc = report("{results}/{sample}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk/{sample}.variantQC.png")

