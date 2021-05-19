#!/home/isshamie/software/anaconda2/envs/vireo/bin python3

import os
from os.path import join, dirname
import papermill as pm
from src.config import ROOT_DIR
import copy
#join(ROOT_DIR, src)



# Flatten dictionary parameters
# if "mttrace" in config:
#     for c in config['mttrace']:
#         config[c] = copy.deepcopy(config['mttrace'][c])
# # Flatten dictionary parameters
# if "mgatk" in config:
#     for c in config['mgatk']:
#         config[c] = copy.deepcopy(config['mgatk'][c])
if "multiplex" in config:
    for c in config['multiplex']:
        if c in config:
            print(c, "Duplicate values in config.")
        config[c] = copy.deepcopy(config['multiplex'][c])


num_cells = config["pseudo_multiplex"]["num_cells"]
is_prop = config["pseudo_multiplex"]["is_proportional"]

# wildcard_constraints:
#     minAF= "[0-9]+"

mt_minC = config["mt"]["minCOUNT"]
mt_minAF = config["mt"]["minMAF"]
minC = config["pre_variant"]["minCOUNT"]
minAF = config["pre_variant"]["minMAF"]
print('minC', minC)

if "prefix" not in config:
    config["prefix"] = ""

rule all:
    input:
           #expand("data/{prefix}/vcf/{name}_cellSNP_minC{minC}_minAF{minAF}", prefix=config["prefix"], name=config["samples"], minC=minC, minAF=minAF),
           #expand("data/{prefix}/vcf/pseudo/minC{minC}_minAF{minAF}/numC{num_cells}_isprop{is_prop}/cellSNP.tag.AD.mtx", prefix=config["prefix"],
           #       num_cells=num_cells, is_prop=is_prop, minC=minC, minAF=minAF),
           #expand("data/{prefix}/vcf/pseudo/minC{minC}_minAF{minAF}/numC{num_cells}_isprop{is_prop}/pseudo.ipynb", prefix=config["prefix"],  name=config["samples"],
           #       num_cells=num_cells, is_prop=is_prop, minC=minC, minAF=minAF),
           #MT
           expand("data/{prefix}/chrM/{name}_cellSNP_minC{mt_minC}_minAF{mt_minAF}", prefix=config["prefix"], name=config["samples"],
                  mt_minC=mt_minC,mt_minAF=mt_minAF),
           #expand("results/{prefix}/chrM/{name}_cellSNP_minC{mt_minC}_minAF{mt_minAF}/lineage_chrM.ipynb", prefix=config["prefix"], name=config["samples"], mt_minC=mt_minC,mt_minAF=mt_minAF),
           expand("data/{prefix}/chrM/output/{name}_cellSNP_minC{mt_minC}_minAF{mt_minAF}/{name}_multiplex.ipynb",prefix=config["prefix"], name=config["samples"], mt_minC=mt_minC,mt_minAF=mt_minAF, num_cells=num_cells),
           expand("data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/multiplex.ipynb",
                  prefix=config["prefix"], name=config["samples"], mt_minC=mt_minC,mt_minAF=mt_minAF, num_cells=num_cells, is_prop=is_prop),
           expand("data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/lineages/lineage.ipynb", prefix=config["prefix"], name=config["samples"],
                  mt_minC=mt_minC,mt_minAF=mt_minAF, num_cells=num_cells, is_prop=is_prop),
           expand("data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/dendrograms/af_dendro.ipynb",
                  prefix=config["prefix"], name=config["samples"], mt_minC=mt_minC,mt_minAF=mt_minAF, num_cells=num_cells, is_prop=is_prop),
           expand("data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/enrichment/.status",
                  prefix=config["prefix"], name=config["samples"], mt_minC=mt_minC,mt_minAF=mt_minAF, num_cells=num_cells, is_prop=is_prop)



def get_bam(wildcards):
    print('bam', os.path.join(config["indir"], wildcards.name, "outs/possorted_bam.bam"))
    return os.path.join(config["indir"], wildcards.name, "outs/possorted_bam.bam")


all_bams = list(map(lambda x: os.path.join(config["indir"], x, "outs/possorted_bam.bam"),
                    config["samples"]))


def get_barcodes(wildcards):
    #"/data2/isshamie/mito_lineage/data/{prefix}/processed/mtscATAC/2020_11_18_Croker/MTblacklist_mtasnucl/PBMC_J/outs/filtered_peak_bc_matrix/barcodes.tsv"
    return os.path.join(config["indir"], wildcards.name,"outs/filtered_peak_bc_matrix/barcodes.tsv")


def get_mt_bam(wildcards):
    if "use_mt_bam" in config and config["use_mt_bam"]:
        return os.path.join(config["mttrace_dir"], wildcards.name, "MT", f"{wildcards.name}.MT.bam")
    else:
        return get_bam(wildcards)



### Option 1: Run cellSNP for chrM, then run vireo, for each
###           sample individually.

### Option 2: Run cellSNP for variants, then combine patients only IF
###           it's pseudobulk. Then run vireo donor notebook on
###           aggregated samples.
########################################################################
## 1. Pseudo multiplex create
########################################################################
rule VCF_cellSNP:
    """Converts 10x output to a cellSNP file. 
    
    Takes in coverage and number of cells thresholds to save memory.
    """

    input: get_bam
    output: directory("data/{prefix}/vcf/{name}_cellSNP_minC{minC}_minAF{minAF}")
    params:
        variants_f=config["variants_f"],
        barcode_f=get_barcodes,
        minMAF=lambda wildcards: wildcards.minAF,
        minCOUNT=lambda wildcards: wildcards.minC,
        exclFLAG=772
    log: "logs/{prefix}/cellSNP/{name}_cellSNP_minC{minC}_minAF{minAF}.log"
    shell:
        "cellsnp-lite -s {input} --exclFLAG {params.exclFLAG} -p 10 --minCOUNT {params.minCOUNT} --minMAF {params.minMAF} -O {output} -R {params.variants_f} -b {params.barcode_f} --printSkipSNPs --UMItag None --cellTAG CB 2> {log} "

rule VCF_pseudo_bulk_create:
    """Subsets and combines the samples to see if it can properly 
        deconvolve the two."""
    #TODO
    #input: lambda wildcards: expand("data/{prefix}/vcf/{name}_cellSNP_{{minC}}_{{minAF}}", name=config["samples"])#, minC=wildcards.minC, minAF=wildcards.minAF)
    input: lambda wildcards: expand("data/{{prefix}}/{name}_cellSNP_minC{{minC}}_minAF{{minAF}}", name=config["samples"])# , minC=wildcards.minC, minAF=wildcards.minAF)
    params:
        num_cells=num_cells,
        is_prop = is_prop,
        outdir = lambda wildcards, output: dirname(output[0])
    output: "data/{prefix}/vcf/pseudo/minC{minC}_minAF{minAF}/numC{num_cells}_isprop{is_prop}/cellSNP.tag.AD.mtx"
    log: "logs/{prefix}/vireo/minC{minC}_minAF{minAF}/numC{num_cells}_isprop{is_prop}.log"
    shell: "python -m src.pseudo_batch {params.outdir} {input} --num_cells {params.num_cells} --is_prop {params.is_prop} > {log} 2>&1"


rule VCF_multiplex_psuedo:
    input:
        AD_F="data/{prefix}/vcf/pseudo/minC{minC}_minAF{minAF}/numC{num_cells}_isprop{is_prop}/cellSNP.tag.AD.mtx", #"data/{prefix}/vireo/pseudo/numC{num_cells}_isprop{is_prop}/cellSNP.tag.AD.mtx",
        notebook=join(ROOT_DIR, "src", "vireo", "vireoSNP_donors.ipynb")
    output: "data/{prefix}/vcf/pseudo/minC{minC}_minAF{minAF}/numC{num_cells}_isprop{is_prop}/pseudo.ipynb"
    params:
        INDIR = lambda wildcards, input: dirname(input.AD_F),
        #DP_F=lambda wildcards, input: join(dirname(input.AD_F), "cellSNP.tag.DP.mtx"),
        #VCF_F=lambda wildcards, input: join(dirname(input.AD_F), "cellSNP.base.vcf.gz")
    shell: "papermill -p INDIR {params.INDIR} {input.notebook} {output}"
########################################################################


########################################################################
## 2. chrMT pileups, lineage-tracing, and pseudo multiplex
########################################################################
rule chrM_cellSNP:
    """Mode 2. Pileup across chrM

    Takes in coverage and number of cells thresholds to save memory.
    """
    input:
        bam_f= get_mt_bam, #get_bam,
        barcode_f=get_barcodes
    output: directory("data/{prefix}/chrM/{name}_cellSNP_minC{mt_minC}_minAF{mt_minAF}")
    log: "logs/{prefix}/cellSNP_chrM/{prefix}/{name}_cellSNP_minC{mt_minC}_minAF{mt_minAF}_chrM.log"
    params:
        minMAF=lambda wildcards: wildcards.mt_minAF,
        minCOUNT=lambda wildcards: wildcards.mt_minC,
        exclFLAG=772 # Read unmapped, not primary alignment, read fails platform qc. Default minMapq is 20
    shell: "cellsnp-lite -s {input.bam_f} --exclFLAG {params.exclFLAG}  -b {input.barcode_f} -O {output} -p 20 --minMAF {params.minMAF} --minCOUNT {params.minCOUNT} --gzip --chrom chrM --printSkipSNPs --UMItag None --cellTAG CB 2> {log}"#chrM"

rule chrM_multiplex_single:
    """Run multiplex not with pseudo-aggregrate, just for one file"""
    input:
        INDIR= "data/{prefix}/chrM/{name}_cellSNP_minC{mt_minC}_minAF{mt_minAF}"
    output:
        out_note="data/{prefix}/chrM/output/{name}_cellSNP_minC{mt_minC}_minAF{mt_minAF}/{name}_multiplex.ipynb",
    params:
        N_DONORS=config["N_DONORS"],
        notebook=join(ROOT_DIR, "src", "vireo", "1_MT_Donors_multiplex_single.ipynb" ),
        OUTDIR = lambda wildcards, output: dirname(output.out_note),
    shell: "papermill -p INDIR {input.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} {params.notebook} {output}"
    # output: "results/{prefix}/chrM/{name}_cellSNP_minC{mt_minC}_minAF{mt_minAF}/lineage_chrM.ipynb" #"data/{prefix}/vireo/{name}_chrM.ipynb"
    # params:
    #     AD_F=lambda wildcards, input: join(input.dat, "cellSNP.tag.AD.mtx"),
    #     DP_F=lambda wildcards, input: join(input.dat, "cellSNP.tag.DP.mtx"),
    #     #VCF_F=lambda wildcards, input: join(input.dat, "cellSNP.base.vcf.gz")
    # shell: "papermill -p AD_F {params.AD_F} -p DP_F {params.DP_F} {input.notebook} {output}"

## ChrMT Pseudo-multiplex
rule chrM_pseudo_bulk_create:
    """Subsets and combines the samples to see if it can properly
        deconvolve the two."""
    #TODO
    input: lambda wildcards: expand("data/{{prefix}}/chrM/{name}_cellSNP_minC{{mt_minC}}_minAF{{mt_minAF}}", name=config["samples"])# , minC=wildcards.minC, minAF=wildcards.minAF)
    params:
        num_cells=num_cells,
        is_prop = is_prop,
        outdir = lambda wildcards, output: dirname(output[0]),
        prefix= ','.join(config["samples"]) # make it as a list
    output: "data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/cellSNP.tag.AD.mtx"
    log: "logs/{prefix}/vireo/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}.log"
    shell: "python -m src.pseudo_batch {params.outdir} {input} --num_cells {params.num_cells} --is_prop {params.is_prop} --samples {params.prefix} > {log} 2>&1"


rule chrM_multiplex_pseudo:
    input:
        AD_F="data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/cellSNP.tag.AD.mtx",
    output:
        out_note="data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/multiplex.ipynb",
    params:
        N_DONORS=config["N_DONORS"],
        notebook=join(ROOT_DIR, "src", "vireo", "1_MT_Donors_multiplex.ipynb" ),
        sample_names= ','.join(config["samples"]), # make it as a list
        INDIR = lambda wildcards, input: dirname(input.AD_F),
        OUTDIR = lambda wildcards, output: dirname(output.out_note),
        #sample_csv = config["sample_csv"]
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p sample_names {params.sample_names} {params.notebook} {output}"


rule chrM_lineage_pseudo:
    input: "data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/multiplex.ipynb"
    output: "data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/lineages/lineage.ipynb"
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["N_DONORS"],
        notebook=join(ROOT_DIR, "src", "vireo", "2_MT_Lineage_Construct.ipynb"),
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} {params.notebook} {output}"


rule chrM_plotAF:
    input:
        #clone="data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/lineages/lineage.ipynb",
        mult="data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/multiplex.ipynb"
    output: "data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/dendrograms/af_dendro.ipynb"
    params:
        #INDIR = lambda wildcards, input: dirname(input[0]),
        INDIR = lambda wildcards, input: dirname(input.mult),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["N_DONORS"],
        sample_names= ','.join(config["samples"]), # make it as a list
        notebook=join(ROOT_DIR, "src", "vireo", "3_MT_Donors_Dendrogram.ipynb"),
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS}  {params.notebook} {output}"


rule chrM_enrichment:
    input:
        "data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/multiplex.ipynb",
        "data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/lineages/lineage.ipynb"
    output: "data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/enrichment/.status"
    params:
        donors_indir = lambda wildcards, input: dirname(input[0]),
        clones_indir = lambda wildcards, input: dirname(input[1]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["N_DONORS"],
        nclones= config["n_clone_list"],
        script=join(ROOT_DIR, "src", "vireo", "lineage_enrichment.py"),
    shell: "python {params.script} {params.donors_indir} {params.clones_indir} {params.OUTDIR} {params.N_DONORS} {params.nclones}"


# rule chrM_pseudo_multiplex:
#     input:
#         AD_F="data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/cellSNP.tag.AD.mtx",
#         notebook=join(ROOT_DIR, "src", "vireo", "vireoSNP_MT_Donors_proper.ipynb" )
#     output: "data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/pseudo.ipynb"
#     params:
#         INDIR = lambda wildcards, input: dirname(input.AD_F),
#         #DP_F=lambda wildcards, input: join(dirname(input.AD_F), "cellSNP.tag.DP.mtx"),
#         #VCF_F=lambda wildcards, input: join(dirname(input.AD_F), "cellSNP.base.vcf.gz")
#     shell: "papermill -p INDIR {params.INDIR} {input.notebook} {output}"
# ################################################################################

# rule deconvolve_individual:
#     input: "data/{prefix}/{name}_cellSNP_minC{minC}_minAF{minAF}"
#     output: "data/{prefix}/vireo_MT_individual_donors/{name}_minC{minC}_minAF{minAF}.ipynb"
#     shell: "papermill -p AD {input}.cellSNP.tag.AD.mtx -p DP {input}.cellSNP.tag.DP.mtx"
#     "src/vireo/vireoSNP_clones.ipynb {output}"


#
# rule cellSNP_all_chr:
#     """Mode 2. Pileup across all positions
#
#     Takes in coverage and number of cells thresholds to save memory.
#     """
#     input:
#         bam_f=get_bam,
#         barcode_f=get_barcodes
#     output: directory("data/{prefix}/{name}_minC{minCOUNT}_minMAF{minMAF}_cellSNP")
#     params:
#           minMAF=lambda wildcards: wildcards.minMAF, #0.01
#           minCOUNT=lambda wildcards: wildcards.minCOUNT #100
#     shell: "cellsnp-lite -s {input.bam_f} -b {input.barcode_f} -O cellsnp_output -p 20 --minMAF {params.minMAF} --minCOUNT {params.minCOUNT} --gzip"


# ########################################
# rule pseudo_bam_create:
#     """Are the cell barcodes the same numbers"""
#     input: all_bams
#     output: "data/{prefix}/vcf/pseudo/mixed.bam"
#     run:
#         count = 1
#         tmp = []
#         to_divide = 1/len(input)*100 #Percent to subsample from each sample
#         to_divide = f"{to_divide:.0f}"
#         for b in input:
#             cmd = f"samtools view -bs 42.{to_divide} {b} > {output}.{count}"
#             count += 1
#             shell(f"{cmd}") # Subsamples. 42 is random int. .1 is 10 percent)
#             tmp.append(f"{output}.{count}")
#         shell(f"samtools merge {output} {' '.join(tmp)}")
#         shell(f"rm {' '.join(tmp)}")
#
#     #shell: "samtools view -bs 42.1 {input} > {output}" # Subsamples. 42 is random int. .1 is 10 percent
# rule pseudo_bulk_cellSNP:
#     input: "data/{prefix}/vcf/pseudo/mixed.bam"
#     output: "data/{prefix}/vcf/pseudo/mixed"
#     params:
#         variants_f=config["variants_f"],
#         barcode_f=get_barcodes
#     log: "logs/{prefix}/psuedobulk_cellsnp.log"
#     shell:
#         "cellsnp-lite -s {input} -p 10 --minCOUNT 20 --minMAF 0.01 -O {output} -R {params.variants_f} -b {params.barcode_f} --printSkipSNPs --UMItag None --cellTAG CB 2> {log}"
# ########################################


# rule pseudo_bulk_run_vireo:
#     """
#     """
#     shell: ""
# rule pseudo_bulk_qc:
#     """
#     """
#     shell: ""
#     #input:
