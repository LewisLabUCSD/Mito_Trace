#from src.config import ROOT_DIR
#workdir: ROOT_DIR
wildcard_constraints:
    cellr='True|False'

import os
from os.path import join, dirname
import pandas as pd
from snakemake.utils import validate
from Bio import SeqIO
import pickle
import subprocess as subp
from datetime import datetime
import copy
import numpy as np
import snakemake

from snakemake.utils import min_version
min_version("6.0")

# Flatten dictionary parameters
if "mttrace" in config:
    for c in config['mttrace']:
        config[c] = copy.deepcopy(config['mttrace'][c])
# Flatten dictionary parameters
if "mgatk" in config:
    for c in config['mgatk']:
        config[c] = copy.deepcopy(config['mgatk'][c])

print('config')
print(config)
cellr_bc = config["mttrace"]["use_cellr_barcode"]
num_reads_filter = config["num_reads_filter"]
maxBP = config["maxBP"]
ref_fa = config["ref_fa"]
#print(pd.read_table(config["samples"], dtype=str,sep=','))
samples = pd.read_table(config["samples"], dtype=str,sep=',').set_index(["sample_name"], drop=False)

ncells_thresh_mgatk=config['ncells_thresh_mgatk']

num_cells = config['multiplex']["pseudo_multiplex"]["num_cells"]
is_prop = config['multiplex']["pseudo_multiplex"]["is_proportional"]

#print('index',samples.index)
res = config["results"]
#mq = config["mapq"]

ft = config["filters"]
#workdir: config["work_dir"]
rule all:
    input:
        expand("{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz",
               results=res,sample=samples["sample_name"].values, num_read=config["num_reads_filter"]),
        expand("{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/mgatk/{sample}.variant.rds", #lowC{low_cov_thresh}_cellT{ncells_thresh_mgatk}.variant.rds",
               results=res,sample=samples["sample_name"].values,
               cellr_bc=cellr_bc, num_read=num_reads_filter, low_cov_thresh=config["low_cov_thresh"], ncells_thresh_mgatk=ncells_thresh_mgatk),
        expand("{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}_MT_position_coverage.png",
               results=res,sample=samples["sample_name"].values, num_read=config["num_reads_filter"], cellr_bc=config["use_cellr_barcode"]),
        expand("{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.variant.rds", #lowC{low_cov_thresh}_cellT{ncells_thresh_mgatk}.variant.rds",
               results=res,sample=samples["sample_name"].values,
               cellr_bc=cellr_bc, num_read=num_reads_filter,
               min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh']), #, low_cov_thresh=config["low_cov_thresh"], ncells_thresh_mgatk=ncells_thresh_mgatk)
        expand("{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx",
               results=res,cellr_bc=cellr_bc, num_read=num_reads_filter,
               min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh']),
        expand("{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/clones.ipynb",
               results=res,cellr_bc=cellr_bc, num_read=num_reads_filter,
               min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh']),
        expand("{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/dendrograms/af_dendro.ipynb",
               results=res,cellr_bc=cellr_bc, num_read=num_reads_filter,
               min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh']),
        expand("{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/enrichment/.status",
               results=res, cellr_bc=cellr_bc, num_read=num_reads_filter,
               min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh']),
        expand("{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb",
               results=res, cellr_bc=cellr_bc, num_read=num_reads_filter,
               min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], n_clones=config['multiplex']['n_clone_list']),
        # expand("{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones_knn/donor{d}_clones.ipynb",
        #        results=res, cellr_bc=cellr_bc, num_read=num_reads_filter,
        #        min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
        #        het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], d=np.arange(config["multiplex"]["N_DONORS"]))
        # # expand("{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/cells_BC.csv",
        #        results=res, cellr_bc=cellr_bc, num_read=num_reads_filter,
        #        min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
        #        het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], n_clones=config['multiplex']['n_clone_list']),
        # expand("{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/clones/clones.ipynb",
        #        results=res,sample=samples["sample_name"].values,
        #        cellr_bc=cellr_bc, num_read=num_reads_filter,
        #        min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
        #        het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh']),



# def get_sample_bai(wildcards):
#     bai = os.path.join(samples.loc[wildcards.sample, "raw"],samples.loc[wildcards.sample, "bam"]) + ".bam.bai"
#     return bai


def get_sample_bam(wildcards):
    bam = samples.loc[wildcards.sample,"bam_f"]
    print(bam)
    return bam # + ".bam"


def get_sample_barcodes(wildcards):
    return samples.loc[wildcards.sample, "barcode_f"]
    #return os.path.join(dirname(samples.loc[wildcards.sample,"bam_f"]), "filtered_feature_bc_matrix/barcodes.tsv")


def results_dir():
    return config["results"]

# rule orig_index:
#     input:  get_raw_bam #"{raw_f}"  #lambda wildcards: f"{config['bam'][wildcards.sample]}"
#     output: "{bam_f}.bam.bai" #"{raw}/{bam}.bai"
#     shell: "samtools index {input}"

# rule cp_bam:
#     """Move bam file to current location"""
#     input: get_sample_bam
#     output:
#         bam = temp("{results}/data/{sample}/00_bam/{sample}.bam"),
#         #bai = "{results}/data/{sample}/00_bam/{sample}.bam.bai"
#     shell: "cp {input} {output}"
# #    run:
#         # shell("ln -s {input} {output.bam}"),
#         # shell("ln -s {input}.bai {output.bai}")
#
# rule index_bam:
#     """Index the bam file"""
#     input: "{results}/data/{sample}/00_bam/{sample}.bam"
#     output: temp("{results}/data/{sample}/00_bam/{sample}.bam.bai")
#     shell: "samtools index {input}"

rule link_bam:
    input: get_sample_bam
    output: "{results}/data/{sample}/00_bam/{sample}.bam"
    shell: 'ln -sr {input} {output}'

rule index_bam:
    """Index the bam file"""
    input: "{results}/data/{sample}/00_bam/{sample}.bam"#rules.link_bam.output
    output: "{results}/data/{sample}/00_bam/{sample}.bam.bai"
    shell: "samtools index {input}"


# rule filter_bam:
#     """Filter reads with low MAPQ and other potential filters """
#     input:
#         get_sample_bam,
#         #in_bam = "{results}/data/{sample}/00_bam/{sample}.bam",
#         #in_bai = "{results}/data/{sample}/00_bam/{sample}.bam.bai",
#     output:
#         mq_bam = "{results}/data/{sample}/01_bam_filter/{sample}.bam",
#         mq_bai = "{results}/data/{sample}/01_bam_filter/{sample}.bam.bai"
#
#     run:
#             if "rm_duplicates" in config and config["rm_duplicates"]:
#                 shell("samtools view {input.in_bam} -F 1284 -q {wildcards.mapq} -h -b > {output.mq_bam}")
#             else:
#                 shell("samtools view {input.in_bam} -q {wildcards.mapq} -h -b > {output.mq_bam}")
#             shell("samtools index {output.mq_bam}")


rule MT_map:
    """Extract the MT genome"""
    # input:
    #     bam = "{results}/data/{sample}/01_bam_filter/{sample}.bam",
    #     bai = "{results}/data/{sample}/01_bam_filter/{sample}.bam.bai"
    input:
        #bam = "{results}/data/{sample}/00_bam/{sample}.bam",
        bai = "{results}/data/{sample}/00_bam/{sample}.bam.bai",
    output:
        mt_bam="{results}/data/{sample}/MT/{sample}.MT.bam",
        mt_bai="{results}/data/{sample}/MT/{sample}.MT.bam.bai"
    params:
        bam = lambda wildcards, input: input.bai.replace('.bai', ''),
        mt_chr=config["mito_character"],

    run:
        #shell("samtools {input.bam}")
        shell("samtools view -b {params.bam} {params.mt_chr} > {output.mt_bam}")
        shell("samtools index {output.mt_bam}")

rule get_bigwig:
    """Extract the MT genome"""
    input: "{results}/data/{sample}/MT/{sample}.MT.bam"
    output:
        coverage="{results}/data/{sample}/MT/{sample}.MT.bw"
    shell: "bamCoverage -b {input} -o {output}"


rule barcode_data:
    """Loop through the bam file and extract the barcode information."""
    input:
        mt_bam="{results}/data/{sample}/MT/{sample}.MT.bam",
        mt_bai="{results}/data/{sample}/MT/{sample}.MT.bam.bai"
    output: "{results}/data/{sample}/MT/{sample}_barcode_data.p",
    params: mt_chr=config["mito_character"]
    shell:
         "python src/bam_barcodes_function.py {input.mt_bam} {output} {params.mt_chr}"


rule barcode_filter:
    input:
        barcode_p = "{results}/data/{sample}/MT/{sample}_barcode_data.p",
        cellr_f = get_sample_barcodes
    output: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_barcode_data.p"
    params: cellr_bc =  lambda wildcards: wildcards.cellr_bc
    shell: "python src/filter_barcodes.py {input} {params} {output}"


## TODO
# rule compare_CB_with_cellranger_CB_list:
#     input:
#         cb_list = "{results}/data/{sample}/MT/{sample}_barcode_data.p",
#         cellranger_list = "{results}/data/{sample}/MT/{sample}.barcodes.txt"
#     output:
#         "{results}/data/{sample}/MT/barcodes_compare/barcodes_detected.png"
#     shell:
#

# rule plot_CB_coverage:
#     """Plot the MT coverage of the single cells"""
#     input: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_barcode_data.p"
#     output: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_CB_coverage_hist.png"
#     shell:
#         "python src/plot_CB_coverage.py {input} {output}"

rule sortCB:
    input:
        mt_bam="{results}/data/{sample}/MT/{sample}.MT.bam",
        mt_bai="{results}/data/{sample}/MT/{sample}.MT.bam.bai"
    output: "{results}/data/{sample}/MT/{sample}.MT.CB.bam"
    shell: "samtools sort -t CB {input.mt_bam} > {output}"


rule scBam:
    """Extract each single-cell and put into respective bam file"""
    input: "{results}/data/{sample}/MT/{sample}.MT.CB.bam"
        #mt_bam="{results}/data/{sample}/MT/{sample}.MT.bam",
        #mt_bai="{results}/data/{sample}/MT/{sample}.MT.bam.bai"
    output: directory("{results}/data/{sample}/MT/{sample}_scBam")
    threads: 18
    #directory("/data/isshamie/miito_lineage/mttrace/{sample}/MT/{sample}_scBam")
    shell:
        "python src/split_by_CB.py {input} {output}"


rule scPileup:
    """Run the first part of the MT-genotype function by getting the read pileups for each bam file for each nucleotide and overall coverage"""
    input:
        scBam = "{results}/data/{sample}/MT/{sample}_scBam",
        #scBam = "{results}/data/{sample}/MT/{sample}_scBam",
        barcodes = "{results}/data/{sample}/MT/{sample}_barcode_data.p",
    output:
          directory("{results}/data/{sample}/MT/{sample}_scPileup_{num_read}")
    params:
        base_quality = config['base_quality']
    threads: 18
    shell:
         "python src/scPileup_counts.py {input.scBam} {output} {input.barcodes} {wildcards.num_read} {params.base_quality}"


def concat_files(directory, samplename, nt):
    cmd = f"find {directory} -type f -name *.{nt}.txt -exec cat {{}} \; > {samplename}_all.{nt}.txt"
    print(cmd)
    subp.check_call(str(cmd), shell=True)
    cmd = f"find {directory} -type f -name *.{nt}.minus.txt -exec cat {{}} \; > {samplename}_all.{nt}.minus.txt"
    print(cmd)
    subp.check_call(str(cmd), shell=True)
    return


rule scConcat:
    input:
        scPileup_dir = "{results}/data/{sample}/MT/{sample}_scPileup_{num_read}",
        scBam = "{results}/data/{sample}/MT/{sample}_scBam"
    output:
        "{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.minus.txt"
    threads: 28
    params:
        samplename = lambda wildcards, output: output[0].split("_all.coverage.minus.txt")[0]
          #"data/processed/{sample}/scPileup_concat/{sample}_{num_read}"
    run:
        for n in ["A", "C", "G", "T", "coverage"]:
            concat_files(input.scPileup_dir, params.samplename, n)


rule scPileup_concat_strands:
    """ Run the second part of the MT-genotype pipeline, which just concatenates all the pileup data for each nucleotide and overall."""
    input:  "{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.minus.txt"
    output:
        all = "{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz"
    params:
        concat_dir = lambda wildcards, input: dirname(input[0]),
        samplename = lambda wildcards, output: output.all.split("_all.coverage.strands.txt.gz")[0]
          #"data/processed/{sample}/scPileup_concat/{sample}_{num_read}"
    shell:
         "python src/scPileup_concat.py {params.concat_dir} {params.samplename}"

rule plot_CB_coverage:
    """Plot the MT coverage of the single cells"""
    input:
         all = "{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz",
         barcodes = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_barcode_data.p"
    output: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_CB_coverage_hist_minReads{num_read}.png"
    shell:
        "python src/plot_CB_coverage.py {input} {output}"


rule scPileup_MT_matrix:
    """Create the position-by-cell coverage matrix"""
    input:
        all = "{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz",
        barcode_p = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_barcode_data.p"
    output:
        sc_coverage_f = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/sc_coverage.csv"
    threads: 16
    shell:
        "python src/plot_heatmap_coverage.py sc_mt {input.barcode_p} {input.all} {output.sc_coverage_f} {maxBP}"


def run_filter_cell(in_f, prefix, barcode_p, n=""):
    df = pd.read_csv(in_f)
    df["CB"] = df["CB"].str.replace(".bam","")
    barcodes = pickle.load(open(barcode_p, "rb"))
    if isinstance(barcodes, dict):
        df = df[df["CB"].isin(list(barcodes.keys()))]
    else:
        df = df[df["CB"].isin(barcodes)]
    if n == "coverage":
        df = df.iloc[:,:3]
    df.to_csv(prefix+f".{n}.strands.txt.gz", header=None, index=None, compression='gzip')
    return


rule filter_cell_bc:
    """Extracts only the relevant cell barcodes and removes the .bam from the barcode names."""
    input:
        all = "{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz",
        barcode_p = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_barcode_data.p"
    output:
        "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.coverage.strands.txt.gz"
    run:
        for n in ["A", "C", "G", "T", "coverage"]:
            curr_f = input.all
            curr_f = curr_f.replace(".coverage.", "." + n + ".")
            df = pd.read_csv(curr_f)
            df["CB"] = df["CB"].str.replace(".bam","")
            barcodes = pickle.load(open(input.barcode_p, "rb"))
            if isinstance(barcodes, dict):
                df = df[df["CB"].isin(list(barcodes.keys()))]
            else:
                df = df[df["CB"].isin(barcodes)]
            curr_out_f = output[0]
            curr_out_f = curr_out_f.replace(".coverage.", "." + n + ".")
            if n == "coverage":
                df = df.iloc[:,:3]
            df = df.sort_values(["CB", "Position"])
            df.to_csv(curr_out_f, header=None, index=None, compression='gzip')

def get_ref(wildcards):
    w = wildcards
    mito = config['mito_character']
    return f"{w.results}/{w.sample}/mapq_{w.mapq}/cellr_{w.cellr_bc}/{w.sample}_{w.num_read}/{mito}_refAllele.txt"


rule plot_scPileup_MT_matrix:
    """Plot the posiitonal coverages."""
    input:
        sc_coverage_f = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/sc_coverage.csv"
    output:
        save_f_heat = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}_MT_position.png",
        save_f_coverage = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}_MT_position_coverage.png"
    shell:
        #"python src/plot_heatmap_coverage.py plot {input.sc_coverage_f} {output.save_f_coverage}"
        "python src/plot_heatmap_coverage.py plot {input.sc_coverage_f} {output.save_f_heat} {output.save_f_coverage}"


def get_filt(w):
    return w.min_cells, w.min_reads, w.topN, w.het_thresh, w.min_het_cells, w.het_count_thresh, w.bq_thresh


rule create_filters:
    input:
        concat_dir = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.coverage.strands.txt.gz"
    #output:  "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/af_by_cell.tsv"
    output:
        af_f = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/af_by_cell.tsv",
        cov = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt"
    params:
        concat_d = lambda wildcards, input: dirname(input.concat_dir),
        ref_fa = config['mt_ref_fa'],
        name = lambda wildcards: wildcards.sample,
        filt_params = get_filt,
    resources:
        mem_mb=90000
    #log: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/.outcfg"
    shell: "python src/calculate_AF_by_cell.py {params.concat_d} {output.af_f} {params.ref_fa} {params.name} {params.filt_params}"# --log {log}"


rule get_refAllele:
    #input: config["mt_ref_fa"],
    params: config["chrM_refAllele"]
    output: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/chrM_refAllele.txt"
    shell: 'cp {params} {output}'


rule mgatk:
    """ Run both toSeurat and call variants in one script"""
    input:
        all = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt",
        refAllele = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/chrM_refAllele.txt"
    output:
        vars_f = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.variant.rds",
        vars_qc = report("{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.variantQC.png")
    params:
        data_dir=lambda wildcards, input: dirname(input.all),
        sample = lambda wildcards: wildcards.sample,
    shell:
        "./R_scripts/wrap_mgatk.R {params.data_dir} filter_mgatk/{params.sample} FALSE"


rule mgatk_to_vireoIn:
    input: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.variant.rds"
    output: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx"
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        sample = lambda wildcards: wildcards.sample
    shell:
        "python src/mgatk_to_vireo.py {params.indir} {params.outdir} {params.sample}"


rule merged:
    input:
        lambda wildcards: expand("{{results}}/data/{sample}/MT/cellr_{{cellr_bc}}/{sample}_{{num_read}}/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx",
                                        sample=samples["sample_name"].values)# , minC=wildcards.minC, minAF=wildcards.minAF)
    #input: lambda wildcards: expand("data/{{prefix}}/chrM/{name}_cellSNP_minC{{mt_minC}}_minAF{{mt_minAF}}", name=config["samples"])
    params:
        num_cells=num_cells,
        is_prop = is_prop,
        indirs = lambda wildcards, input: [dirname(x) for x in input],
        outdir = lambda wildcards, output: dirname(output[0]),
        prefix= ','.join(samples["sample_name"])
    #log: "logs/{prefix}/vireo/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}.log"
    output: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx"
    shell: "python -m src.pseudo_batch {params.outdir} {params.indirs} --num_cells {params.num_cells} --is_prop {params.is_prop} --samples {params.prefix}" # > {log} 2>&1"

# from snakemake.utils import min_version
# min_version("6.0")

# module other_workflow:
#     snakefile: "other_workflow/Snakefile"
#
# use rule * from other_workflow as other_*
rule multiplex:
    input: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx"#rules.merged.output
    output:
        "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/multiplex.ipynb"
    params:
        N_DONORS=config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "1_MT_Donors_multiplex.ipynb" ),
        sample_names= ','.join(samples["sample_name"].values), # make it as a list
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        workdir = os.getcwd(),
        to_elbo = False
    shell: "papermill --cwd {params.workdir} -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p sample_names {params.sample_names} -p to_elbo {params.to_elbo} {params.notebook} {output}"



# Extract pileups for each donor from original filter matrices
rule scPileup_filter_mgatk:
    input:
        cov = expand("{{results}}/data/{sample}/MT/cellr_{{cellr_bc}}/{sample}_{{num_read}}/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/{sample}.coverage.txt", sample=samples['sample_name'].values),
        mult = "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/multiplex.ipynb"
    output:
        cov = expand("{{results}}/data/merged/MT/cellr_{{cellr_bc}}/numread_{{num_read}}/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/d{d}.coverage.strands.txt.gz",
                     d=config["multiplex"]["N_DONORS"])
    params:
        indir = lambda wildcards, input: dirname(input.cov),
        outdir = lambda wildcards, output: dirname(output.cov),
        cells_meta = lambda wildcards, input: join(dirname(input.mult), "cells_meta.csv"),
        samples = lambda wildcards: wildcards.sample

    run:
        cells_meta = pd.read_csv(params.cells_meta, sep='\t')
        for d, df in cells_meta.groupby('donor'):
            nt_pileups = {}
            cond_positions = {}
            nt_pileups["coverage"] = pd.DataFrame()

            for ind, curr_cov_f in enumerate(input.cov):
                curr_samp_cov = pd.read_csv(curr_cov_f, header=None) #load sample coverage
                condition = params.samples[ind]
                curr_samp_donor_df = df[df["condition"]==condition]
                filt_df = curr_samp_donor_df[curr_samp_cov[1].isin(curr_samp_donor_df["raw ID"])]
                filt_df[1] = filt_df[1] + "_" + params.sample[ind]
                cond_positions[params.sample[ind]] = set(filt_df[0].values)
                nt_pileups["coverage"] = nt_pileups["coverage"].append(filt_df, ignore_index=True)
            # Filter by overlapping positions:
            pos_to_keep = list(cond_positions.values())[0].intersection(*list(cond_positions.values()))
            nt_pileups["coverage"] = nt_pileups["coverage"][nt_pileups["coverage"][0].isin(pos_to_keep)]

            # Save coverage and depth
            nt_pileups["coverage"].to_csv(output.cov[d])
            depth = nt_pileups["coverage"].groupby(1).sum()[2]/16569
            depth.to_csv(join(dirname(output.cov[d]), f"d{d}.depthTable.txt"), sep='\t', header=False)

            # Filter for the nucleotides as well
            for nuc in ["C", "A", "G", "T"]:
                nt_pileups[nuc] = pd.DataFrame()
                for ind, curr_cov_f in enumerate(input.cov):
                    curr_samp_cov = pd.read_csv(join(dirname(curr_cov_f), f"{params.samples[ind]}.{nuc}.txt"), header=None) #load sample coverage
                    curr_samp_donor_df = df[df["condition"]==condition]
                    filt_df = curr_samp_donor_df[curr_samp_cov[1].isin(curr_samp_donor_df["raw ID"])] # Filter by donor
                    filt_df[1] = filt_df[1] + "_" + params.sample[ind]
                    filt_df = filt_df[filt_df[1].isin(curr_samp_donor_df["ID"])] # Filter by cell ID from multiplex cells_meta
                    filt_df = filt_df[filt_df[0].isin(pos_to_keep)]
                    nt_pileups[nuc] = nt_pileups[nuc].append(filt_df, ignore_index=True)
                # Save coverage
                nt_pileups[nuc].to_csv(join(dirname(output.cov[d]), f"d{d}.{nuc}.txt"))
            #shell: "python donor_filter_mgatk.py {params.indir} {params.cells_meta} {params.outdir}"


rule donor_get_refAllele:
    #input: config["mt_ref_fa"],
    params: config["chrM_refAllele"]
    output: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/chrM_refAllele.txt"
    shell: 'cp {params} {output}'


rule donor_mgatk:
    input: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/d{d}.coverage.strands.txt.gz"
    output: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/donor_mgatk/d{d}.variant.rds"
    params:
        data_dir =lambda wildcards, input: dirname(input[0]),
        donor = lambda wildcards: f"d{wildcards.d}",
    shell:
        "./R_scripts/wrap_mgatk.R {params.data_dir} donor_mgatk/{params.donor} FALSE"


rule donor_mgatk_to_vireoIn:
    input: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/donor_mgatk/d{d}.variant.rds"
    output: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/donor_mgatk/vireoIn/cellSNP.tag.AD.mtx"
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        donor = lambda wildcards: f"d{wildcards.d}"
    shell:
        "python src/mgatk_to_vireo.py {params.indir} {params.outdir} {params.donor}"


rule clones_donor_mgatk:
    input: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/donor_mgatk/vireoIn/cellSNP.tag.AD.mtx"
    output: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/donors_mgatk_in/donor{d}/donor_mgatk/vireoIn/clones/clones.ipynb"
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "2_MT_Lineage_Construct.ipynb"),
        workdir = os.getcwd(),
        n_clone=",".join([str(x) for x in config['multiplex']['n_clone_list']])
    shell: "papermill --cwd {params.workdir} -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p n_clone_list {params.n_clone} {params.notebook} {output}"


rule clones:
    input: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/multiplex.ipynb",
    output: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/clones.ipynb"
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "2_MT_Lineage_Construct.ipynb"),
        workdir = os.getcwd(),
        n_clone=",".join([str(x) for x in config['multiplex']['n_clone_list']])
    shell: "papermill --cwd {params.workdir} -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p n_clone_list {params.n_clone} {params.notebook} {output}"


# from snakemake.utils import R
# rule merge_mgatk:
#     input:
#         expand("{{results}}/data/{sample}/MT/cellr_{{cellr_bc}}/{sample}_{{num_read}}/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/{sample}.variant.rds",
#                sample=samples["sample_name"].values)
#     output:
#         mgatk="{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/mgatk.variants.merged.rds",
#
#     run:
#         R("""
#           SE_list = lapply(input, readRDS)
#           writeRDS(cbind(SE_list), {output.mgatk})
#         """)
#
# rule clones_knn:
#     input:
#         mult="{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/multiplex.ipynb",
#         mgatk="{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/mgatk.variants.merged.rds"
#     output:
#         "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones_knn/donor{d}_clones.ipynb"
#     params:
#         cells_f = lambda wildcards, input: join(dirname(input.mgatk),"cells_meta.tsv"),
#         outdir = lambda wildcards, output: dirname(output[0]),
#         notebook = join("R_scripts", "annotations/call_clones.ipynb"),
#         cells_col = "donor"
#     shell:
#         "papermill -p mgatk_in {input.mgatk} -p outdir {params.outdir} -p cells_col {params.cells_col} {params.notebook} {output}  && jupyter nbconvert --to pdf {output}"


def get_comparisons(wildcards):
    if "comparisons" in config["multiplex"]:
        return f"--tests {config['comparisons']}"
    return ""


rule enrichment:
    input:
        "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/multiplex.ipynb",
        "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/clones.ipynb",
    output: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/enrichment/.status"
    params:
        donors_indir = lambda wildcards, input: dirname(input[0]),
        clones_indir = lambda wildcards, input: dirname(input[1]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS = config["multiplex"]["N_DONORS"],
        nclones = config['multiplex']["n_clone_list"],
        script = join("src",  "lineage_enrichment.py"),
        samples = ",".join(config["multiplex"]['samples']),
        tests = get_comparisons
    shell: "python {params.script} {params.clones_indir} {params.OUTDIR} {params.nclones} {params.samples} {params.tests}"


rule plotAF:
    input:
        "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/multiplex.ipynb"
    output:
        "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/dendrograms/af_dendro.ipynb",
        #report("{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/dendrograms/figures/")
        report(expand("{{results}}/data/merged/MT/cellr_{{cellr_bc}}/numread_{{num_read}}/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/dendrograms/figures/donor{n}_dendrogram.png",
               n=np.arange(config["multiplex"]["N_DONORS"])))
    params:
        #INDIR = lambda wildcards, input: dirname(input[0]),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["multiplex"]["N_DONORS"],
        sample_names= ','.join(samples["sample_name"].values), # make it as a list
        notebook=join("src", "vireo", "3_MT_Donors_Dendrogram.ipynb"),
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} {params.notebook} {output[0]}"


rule plotAF_clones:
    input:
        #clone="data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/lineages/lineage.ipynb",
        mult = "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/multiplex.ipynb"
    output:
        "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/dendrograms/dendrograms.ipynb",
        #report("{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/nfigures/donor{n}_dendrogram.png")
        report(expand("{{results}}/data/merged/MT/cellr_{{cellr_bc}}/numread_{{num_read}}/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/clones/nclones_{nclones}/dendrograms/figures/donor{n}_dendrogram.png",
                n=np.arange(config["multiplex"]["N_DONORS"]),
                nclones=config['multiplex']["n_clone_list"]))
    params:
        #INDIR = lambda wildcards, input: dirname(input[0]),
        INDIR = lambda wildcards, input: dirname(input.mult),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["multiplex"]["N_DONORS"],
        sample_names= ','.join(samples["sample_name"].values), # make it as a list
        notebook=join("src", "vireo", "3_MT_Donors_Dendrogram.ipynb"),
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} {params.notebook} {output[0]}"



rule clones_type_variants:
    input: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/clones.ipynb"#rules.clones.output
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
    output: "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb"#, n_clones=config['multiplex']["n_clone_list"])
    #output: report(lambda wildcards: expand("{{results}}/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb", n_clones=config['multiplex']["n_clone_list"]))
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p n_clones {params.n_clones} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} -p var_thresh {params.var_thresh} -p vars_to_plot {params.vars_to_plot} {params.notebook} {output} && jupyter nbconvert --to pdf {output}"


# rule clones_exp:
#     input: rules.multiplex.output
#     output: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/clones/clones.ipynb"
#     params:
#         INDIR = lambda wildcards, input: dirname(input[0]),
#         OUTDIR = lambda wildcards, output: dirname(output[0]),
#         N_DONORS=config["multiplex"]["N_DONORS"],
#         notebook=join("src", "vireo", "2_MT_Lineage_Construct_separateCond.ipynb"),
#         workdir = os.getcwd(),
#         exp = lambda wildcards: wildcards.sample
#     shell: "papermill --cwd {params.workdir} -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p exp {params.exp} -p N_DONORS {params.N_DONORS} {params.notebook} {output}"


def get_aggr(in_d):
    #bc = join(in_d, 'aggregate', 'outs', 'filtered_peak_bc_matrix', 'barcodes.tsv')
    aggr_csv = join(in_d, 'aggr.csv')
    return aggr_csv


rule merge_lineage_nuclear_barcodes:
    input:
        cl="{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/clones.ipynb",
        aggr=get_aggr(config['mtscATAC_OUTDIR'])
    params:
        N_DONORS=config["multiplex"]["N_DONORS"],
        #n_clones = lambda wildcards: wildcards.n_clones,
        cells_meta = lambda wildcards, input: join(dirname(input[0]),f'lineage{wildcards.n_clones}', "cells_meta.tsv"),
        #n_clones=config['multiplex']['n_clone_list']
    output:
        bc="{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/cells_BC.csv",
        meta="{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/cells_meta_aggr.tsv"
    run:
        cells_meta = pd.read_csv(params.cells_meta, sep='\t')
        aggr = pd.read_csv(input[1])
        #aggr.index = aggr.index + 1 #1-based index
        print('aggr')
        print(aggr)
        samples_d = {val:i+1 for i, val in enumerate(aggr['library_id'].values)}
        cells_meta = cells_meta.fillna('Unassigned')
        cells_meta['BC aggr'] = cells_meta.apply(lambda x: f"{x['raw ID'].replace('-1','')}-{samples_d[x['condition']]}", axis=1)
        cells_meta['lineage ID'] = cells_meta.apply(lambda x: f"{x['donor']}_{x['lineage']}", axis=1)
        cells_meta.to_csv(output[1], sep='\t', index=False)
        cells_meta[['BC aggr', 'lineage ID']].to_csv(output[0], index=False, header=False)
        # Do for each donor
        for d, df in cells_meta.groupby("donor"):
            df[['BC aggr', 'lineage ID']].to_csv(f"{output[0]}_donor{d}.csv", index=False, header=False)

# report(expand("clones/figures/donors_{donor}/nclones_{nclones}/scatterFisherEnrich.png"))
#import snakemake
#snakemake.io.wo
#
#
# rule variants:
#     input: rules.clones.output
#     output:
#         "{results}/data/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/variants/",
#         report(directory("{results}/data/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/variants/figures"))
#     shell: ""
#

# B. Run mgatk before the filters
rule raw_get_depth:
    input:  "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.coverage.strands.txt.gz" #"{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/sc_coverage.csv"
    output: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.depthTable.txt"
    run:
        df = pd.read_csv(input[0], header=None)
        depth = df.groupby(1).sum()[2]/16569
        depth.to_csv(output[0], sep='\t', header=False)


rule raw_get_refAllele:
    input: config["mt_ref_fa"],
    output: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read,[0-9]+}/chrM_refAllele.txt"
    run:
        records = list(SeqIO.parse(input[0], "fasta"))
        mt_seq = records[0].seq
        for ind, val in enumerate(mt_seq):
            if ind == 0:
                ref_str = f"{ind+1}\t{val}"
            else:
                ref_str = f"{ref_str}\n{ind+1}\t{val}"
        with open(output[0], "w") as f:
            f.write(ref_str)


rule raw_mgatk:
    """ Run both toSeurat and call variants in one script"""
    input:
        all = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.coverage.strands.txt.gz",
        depth = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.depthTable.txt",
        refAllele = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read,[0-9]+}/chrM_refAllele.txt"
    output: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/mgatk/{sample}.variant.rds"
    params:
        data_dir=lambda wildcards, input: dirname(input.all),
        sample = lambda wildcards: wildcards.sample,
    shell:
        "./R_scripts/wrap_mgatk.R {params.data_dir} mgatk/{params.sample} TRUE"
