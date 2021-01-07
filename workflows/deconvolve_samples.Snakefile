import os
from os.path import join
import papermill as pm
from src.config import ROOT_DIR


rule all:
    input: expand("data/{name}_cellSNP", name=config["samples"]),
           expand("data/{name}_cellSNP_chrM", name=config["samples"]),
           expand("data/vireo/{name}_chrM.ipynb", name=config["samples"])
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
    output: directory("data/{name}_cellSNP")
    params:
        variants_f=config["variants_f"],
        barcode_f=get_barcodes
    log: "logs/cellSNP/{name}.log"
    shell:
        "cellsnp-lite -s {input} -p 10 --minCOUNT 20 --minMAF 0.01 -O {output} -R {params.variants_f} -b {params.barcode_f} --printSkipSNPs --UMItag None --cellTAG CB 2> {log} "

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
    output: directory("data/{name}_cellSNP_chrM")
    log: "logs/cellSNP_chrMT/{name}.log"
    shell: "cellsnp-lite -s {input.bam_f} -b {input.barcode_f} -O {output} -p 20 --minMAF 0.01 --minCOUNT 100 --gzip --chrom chrM --printSkipSNPs --UMItag None --cellTAG CB 2> {log}"#chrM"


rule vireo_chrM:
    input:
        dat="data/{name}_cellSNP_chrM",
        notebook=join(ROOT_DIR, "src", "vireo", "vireoSNP_donors.ipynb" )
    output: "data/vireo/{name}_chrM.ipynb"
    params:
        AD_F=lambda wildcards, input: join(input.dat, "cellSNP.tag.AD.mtx"),
        DP_F=lambda wildcards, input: join(input.dat, "cellSNP.tag.DP.mtx")
    shell: "papermill -p AD_F {params.AD_F} -p DP_F {params.DP_F} {input.notebook} {output}"


### Option 1: Run cellSNP for chrM, then run vireo, for each
###           sample individually.

### Option 2: Run cellSNP for variants, then combine patients only IF
###           it's pseudobulk. Then run vireo donor notebook on
###           aggregated samples.
rule deconvolve_individual:
    input: "data/{name}_cellSNP"
    output: "data/vireo_donors/{name}.ipynb"
    shell: "papermill -p AD {input}.cellSNP.tag.AD.mtx -p DP {input}.cellSNP.tag.DP.mtx"
    "src/vireo/vireoSNP_clones.ipynb {output}"


rule pseudo_bulk_create:
    """Subsets and combines the samples to see if it can properly 
        deconvolve the two."""
    #TODO
    input: expand("data/{name}_cellSNP", name=config["samples"])
    params: num_reads_total=config["pseudo-multiplex"]["Number Reads Total"]
    output: directory("data/pseudo_bulk/Reads{params.num_reads_total}_subset")
    shell: "python src/pseudo_batch.py {input} {params} {output}"


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



rule pseudo_bulk_run_vireo:
    """
    """
    shell: ""
rule pseudo_bulk_qc:
    """
    """
    shell: ""
    #input:
