#!/home/isshamie/software/anaconda2/envs/vireo/bin python3

import os
from os.path import join, dirname
import papermill as pm
from src.config import ROOT_DIR
#join(ROOT_DIR, src)

num_cells = config["pseudo_multiplex"]["num_cells"]
is_prop = config["pseudo_multiplex"]["is_proportional"]

# wildcard_constraints:
#     minAF= "[0-9]+"

mt_minC = config["mt"]["minCOUNT"]
mt_minAF = config["mt"]["minMAF"]
minC = config["pre_variant"]["minCOUNT"]
minAF = config["pre_variant"]["minMAF"]
print('minC', minC)
rule all:
    input: expand("data/{name}_cellSNP_minC{minC}_minAF{minAF}", name=config["samples"], minC=minC, minAF=minAF),
           expand("data/chrM/{name}_cellSNP_minC{mt_minC}_minAF{mt_minAF}", name=config["samples"],
                  mt_minC=mt_minC,mt_minAF=mt_minAF),
           expand("data/chrM/{name}/cellSNP_minC{mt_minC}_minAF{mt_minAF}/lineage_chrM.ipynb", name=config["samples"], mt_minC=mt_minC,mt_minAF=mt_minAF ),
           expand("data/pseudo/minC{minC}_minAF{minAF}/numC{num_cells}_isprop{is_prop}/cellSNP.tag.AD.mtx",
                  num_cells=num_cells, is_prop=is_prop, minC=minC, minAF=minAF),
           expand("data/pseudo/minC{minC}_minAF{minAF}/numC{num_cells}_isprop{is_prop}/pseudo.ipynb", name=config["samples"],
                  num_cells=num_cells, is_prop=is_prop, minC=minC, minAF=minAF)
           #expand("data/pseudo/mixed.bam")
            #data/pseudo/mixed.bam


def get_bam(wildcards):
    return os.path.join(config["indir"], wildcards.name, "outs/possorted_bam.bam")


all_bams = list(map(lambda x: os.path.join(config["indir"], x, "outs/possorted_bam.bam"),
                    config["samples"]))


def get_barcodes(wildcards):
    #"/data2/isshamie/mito_lineage/data/processed/mtscATAC/2020_11_18_Croker/MTblacklist_mtasnucl/PBMC_J/outs/filtered_peak_bc_matrix/barcodes.tsv"
    return os.path.join(config["indir"], wildcards.name,"outs/filtered_peak_bc_matrix/barcodes.tsv")


def get_mt_bam(wildcards):
    if "use_mt_bam" in config and config["use_mt_bam"]:
        return os.path.join(config["mttrace_dir"], wildcards.name, "mapq_0", f"{wildcards.name}.MT.bam")
    else:
        return get_bam(wildcards)


rule cellSNP:
    """Converts 10x output to a cellSNP file. 
    
    Takes in coverage and number of cells thresholds to save memory.
    """

    input: get_bam
    output: directory("data/{name}_cellSNP_minC{minC}_minAF{minAF}")
    params:
        variants_f=config["variants_f"],
        barcode_f=get_barcodes,
        minMAF=lambda wildcards: wildcards.minAF,
        minCOUNT=lambda wildcards: wildcards.minC,
        exclFLAG=772

    log: "logs/cellSNP/{name}_cellSNP_minC{minC}_minAF{minAF}.log"
    shell:
        "cellsnp-lite -s {input} --exclFLAG {params.exclFLAG} -p 10 --minCOUNT {params.minCOUNT} --minMAF {params.minMAF} -O {output} -R {params.variants_f} -b {params.barcode_f} --printSkipSNPs --UMItag None --cellTAG CB 2> {log} "

#
# rule cellSNP_all_chr:
#     """Mode 2. Pileup across all positions
#
#     Takes in coverage and number of cells thresholds to save memory.
#     """
#     input:
#         bam_f=get_bam,
#         barcode_f=get_barcodes
#     output: directory("data/{name}_minC{minCOUNT}_minMAF{minMAF}_cellSNP")
#     params:
#           minMAF=lambda wildcards: wildcards.minMAF, #0.01
#           minCOUNT=lambda wildcards: wildcards.minCOUNT #100
#     shell: "cellsnp-lite -s {input.bam_f} -b {input.barcode_f} -O cellsnp_output -p 20 --minMAF {params.minMAF} --minCOUNT {params.minCOUNT} --gzip"

rule cellSNP_chrM:
    """Mode 2. Pileup across chrM

    Takes in coverage and number of cells thresholds to save memory.
    """
    input:
        bam_f= get_mt_bam, #get_bam,
        barcode_f=get_barcodes
    output: directory("data/chrM/{name}_cellSNP_minC{mt_minC}_minAF{mt_minAF}")
    log: "logs/cellSNP_chrM/{name}_cellSNP_minC{mt_minC}_minAF{mt_minAF}_chrM.log"
    params:
        minMAF=lambda wildcards: wildcards.mt_minAF,
        minCOUNT=lambda wildcards: wildcards.mt_minC,
        exclFLAG=772
    shell: "cellsnp-lite -s {input.bam_f} --exclFLAG {params.exclFLAG}  -b {input.barcode_f} -O {output} -p 20 --minMAF {params.minMAF} --minCOUNT {params.minCOUNT} --gzip --chrom chrM --printSkipSNPs --UMItag None --cellTAG CB 2> {log}"#chrM"


rule vireo_chrM:
    input:
        dat= "data/chrM/{name}_cellSNP_minC{mt_minC}_minAF{mt_minAF}",
        notebook=join(ROOT_DIR, "src", "vireo", "vireoSNP_clones.ipynb")
    output: "data/chrM/{name}_cellSNP_minC{mt_minC}_minAF{mt_minAF}/lineage_chrM.ipynb" #"data/vireo/{name}_chrM.ipynb"
    params:
        AD_F=lambda wildcards, input: join(input.dat, "cellSNP.tag.AD.mtx"),
        DP_F=lambda wildcards, input: join(input.dat, "cellSNP.tag.DP.mtx"),
        #VCF_F=lambda wildcards, input: join(input.dat, "cellSNP.base.vcf.gz")
    shell: "papermill -p AD_F {params.AD_F} -p DP_F {params.DP_F} {input.notebook} {output}"


### Option 1: Run cellSNP for chrM, then run vireo, for each
###           sample individually.

### Option 2: Run cellSNP for variants, then combine patients only IF
###           it's pseudobulk. Then run vireo donor notebook on
###           aggregated samples.


# rule deconvolve_individual:
#     input: "data/{name}_cellSNP_minC{minC}_minAF{minAF}"
#     output: "data/vireo_MT_individual_donors/{name}_minC{minC}_minAF{minAF}.ipynb"
#     shell: "papermill -p AD {input}.cellSNP.tag.AD.mtx -p DP {input}.cellSNP.tag.DP.mtx"
#     "src/vireo/vireoSNP_clones.ipynb {output}"


rule pseudo_bulk_create:
    """Subsets and combines the samples to see if it can properly 
        deconvolve the two."""
    #TODO
    #input: lambda wildcards: expand("data/{name}_cellSNP_{{minC}}_{{minAF}}", name=config["samples"])#, minC=wildcards.minC, minAF=wildcards.minAF)
    input: lambda wildcards: expand("data/{name}_cellSNP_minC{{minC}}_minAF{{minAF}}", name=config["samples"])# , minC=wildcards.minC, minAF=wildcards.minAF)
    params:
        num_cells=num_cells,
        is_prop = is_prop,
        outdir = lambda wildcards, output: dirname(output[0])
    output: "data/pseudo/minC{minC}_minAF{minAF}/numC{num_cells}_isprop{is_prop}/cellSNP.tag.AD.mtx"
    log: "logs/vireo/minC{minC}_minAF{minAF}/numC{num_cells}_isprop{is_prop}.log"
    shell: "python -m src.pseudo_batch {params.outdir} {input} --num_cells {params.num_cells} --is_prop {params.is_prop} > {log} 2>&1"


rule vireo_psuedo:
    input:
        AD_F="data/pseudo/minC{minC}_minAF{minAF}/numC{num_cells}_isprop{is_prop}/cellSNP.tag.AD.mtx", #"data/vireo/pseudo/numC{num_cells}_isprop{is_prop}/cellSNP.tag.AD.mtx",
        notebook=join(ROOT_DIR, "src", "vireo", "vireoSNP_donors.ipynb" )
    output: "data/pseudo/minC{minC}_minAF{minAF}/numC{num_cells}_isprop{is_prop}/pseudo.ipynb"
    params:
        INDIR = lambda wildcards, input: dirname(input.AD_F),
        #DP_F=lambda wildcards, input: join(dirname(input.AD_F), "cellSNP.tag.DP.mtx"),
        #VCF_F=lambda wildcards, input: join(dirname(input.AD_F), "cellSNP.base.vcf.gz")
    shell: "papermill -p INDIR {params.INDIR} {input.notebook} {output}"



########################################
rule pseudo_bam_create:
    """Are the cell barcodes the same numbers"""
    input: all_bams
    output: "data/pseudo/mixed.bam"
    run:
        count = 1
        tmp = []
        to_divide = 1/len(input)*100 #Percent to subsample from each sample
        to_divide = f"{to_divide:.0f}"
        for b in input:
            cmd = f"samtools view -bs 42.{to_divide} {b} > {output}.{count}"
            count += 1
            shell(f"{cmd}") # Subsamples. 42 is random int. .1 is 10 percent)
            tmp.append(f"{output}.{count}")
        shell(f"samtools merge {output} {' '.join(tmp)}")
        shell(f"rm {' '.join(tmp)}")

    #shell: "samtools view -bs 42.1 {input} > {output}" # Subsamples. 42 is random int. .1 is 10 percent
rule pseudo_bulk_cellSNP:
    input: "data/pseudo/mixed.bam"
    output: "data/pseudo/mixed"
    params:
        variants_f=config["variants_f"],
        barcode_f=get_barcodes
    log: "logs/psuedobulk_cellsnp.log"
    shell:
        "cellsnp-lite -s {input} -p 10 --minCOUNT 20 --minMAF 0.01 -O {output} -R {params.variants_f} -b {params.barcode_f} --printSkipSNPs --UMItag None --cellTAG CB 2> {log}"
########################################


rule pseudo_bulk_run_vireo:
    """
    """
    shell: ""
rule pseudo_bulk_qc:
    """
    """
    shell: ""
    #input:
