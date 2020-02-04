import os
import pandas as pd
from snakemake.utils import validate

#configfile: "parameters/Zhu_Single_Cell_10X_genomics_Human2_002.yaml"

samples = config["samples"]
raw = config["raw"]
bam = config["bam"]
num_reads_filter = config["num_reads_filter"]
maxBP = config["maxBP"]
ref_fa = config["ref_fa"]
samples = pd.read_table(config["samples"], dtype=str).set_index(["sample"], drop=False)

#samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index


#print(raw["A"])

# raw folders, bam files with different name, and sample name (the new name)
#Problem: Indexing the file, and then wanting to switch the name
# Could do:
# A) copy the initial folder by using the input as the raw bam and then the output is the new name. Then index
# B) conf

#

RAW_SAMPLES = samples.apply(lambda x: os.path.join(x["raw"], x["bam"]), axis=1)

print('index',samples.index)

rule all:
    input:
        expand("figures/{sample}/{sample}_CB_coverage_hist.png",sample=samples["sample"].values),
        expand("data/processed/{sample}/{sample}_scPileup_{num_read}",sample=samples["sample"].values, num_read=config["num_reads_filter"]),
        expand("data/processed/{sample}/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.txt.gz",sample=samples["sample"].values, num_read=config["num_reads_filter"]),
        expand("figures/{sample}/{sample}_{num_read}_MT_position.png", sample=samples["sample"].values, num_read=config["num_reads_filter"]),
        expand("figures/{sample}/{sample}_{num_read}_MT_position_coverage.png", sample=samples["sample"].values, num_read=config["num_reads_filter"])

#sample=config["samples"]),
#expand("{raw_f}.bam.bai", raw_f=RAW_SAMPLES),


def get_raw_bai(wildcards):
    print(wildcards.sample)
    bam = os.path.join(samples.loc[wildcards.sample, "raw"],samples.loc[wildcards.sample, "bam"])
    return bam + ".bam.bai"


def get_raw_bam(wildcards):
    print('raw',wildcards.raw)
    print('bam',wildcards.bam)
    bam = os.path.join(samples.loc[wildcards.raw],samples.loc[wildcards.bam]+".bam")
    return bam + ".bam"


def get_sample_bai(wildcards):
    bai = os.path.join(samples.loc[wildcards.sample, "raw"],samples.loc[wildcards.sample, "bam"]) + ".bam.bai"
    return bai


def get_sample_bam(wildcards):
    bam = samples.loc[wildcards.sample,"bam_f"]
    print(bam)
    return bam + ".bam"

# rule orig_index:
#     input:  get_raw_bam #"{raw_f}"  #lambda wildcards: f"{config['bam'][wildcards.sample]}"
#     output: "{bam_f}.bam.bai" #"{raw}/{bam}.bai"
#     shell: "samtools index {input}"


rule cp_bam:
    """Move bam file to current location"""
    input: get_sample_bam
    output: "data/processed/{sample}/{sample}.bam"
    shell: "cp {input} {output}"

rule index_bam:
    """Index the bam file"""
    input: "data/processed/{sample}/{sample}.bam"
    output: "data/processed/{sample}/{sample}.bam.bai"
    shell: "samtools index {input}"

rule MT_map:
    """Extract the MT genome"""
    input:
        bam = "data/processed/{sample}/{sample}.bam",
        bai = "data/processed/{sample}/{sample}.bam.bai"
    output:
        mt_bam="data/processed/{sample}/{sample}.MT.bam",
        mt_bai="data/processed/{sample}/{sample}.MT.bam.bai"
    run:
        #shell("samtools {input.bam}")
        shell("samtools view -b {input.bam} MT > {output.mt_bam}")
        shell("samtools index {output.mt_bam}")


rule barcode_data:
    """Loop through the bam file and extract the barcode information."""
    input:
        mt_bam="data/processed/{sample}/{sample}.MT.bam",
        mt_bai="data/processed/{sample}/{sample}.MT.bam.bai"
    output: "data/processed/{sample}/{sample}_barcode_data.p"
    shell:
         "python src/bam_barcodes_function.py {input.mt_bam} {output}"


rule plot_CB_coverage:
    """Plot the MT coverage of the single cells"""
    input: "data/processed/{sample}/{sample}_barcode_data.p"
    output: "figures/{sample}/{sample}_CB_coverage_hist.png"
    shell:
        "python src/plot_CB_coverage.py {input} {output}"


rule scBam:
    """Extract each single-cell and put into respective bam file"""
    input:
        mt_bam="data/processed/{sample}/{sample}.MT.bam",
        mt_bai="data/processed/{sample}/{sample}.MT.bam.bai"
    output: directory("data/processed/{sample}/{sample}_scBam")
    shell:
        "python src/split_by_CB.py {input.mt_bam} {output}"


rule scPileup:
    """Run the first part of the MT-genotype function by getting the read pileups for each bam file for each nucleotide and overall coverage"""
    input:
        scBam = "data/processed/{sample}/{sample}_scBam",
        barcodes = "data/processed/{sample}/{sample}_barcode_data.p",
    output:
          directory("data/processed/{sample}/{sample}_scPileup_{num_read}")
    # params:
    #     num_reads_filter = config["num_reads_filter"]
    shell:
         "python src/scPileup_counts.py {input.scBam} {output} {input.barcodes} {wildcards.num_read}"


# rule scFilter:
#     input: "data/processed/{sample}/{sample}_scPileup_{num_read}"
#     output:
#         directory("data/processed/{sample}/{sample}_scPileup_Filter_{num_read}")
#     shell: "python src/scFilter.py {input} {wildcards.num_read} {output}"
#

rule scPileup_concat:
    """ Run the second part of the MT-genotype pipeline, which just concatenates all the pileup data for each nucleotide and overall."""
    input:
        scPileup_dir = "data/processed/{sample}/{sample}_scPileup_{num_read}"
    output:
        all = "data/processed/{sample}/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.txt.gz"
    params:
        samplename = lambda wildcards, output: output.all.split("_all.coverage.txt.gz")[0]
          #"data/processed/{sample}/scPileup_concat/{sample}_{num_read}"
    shell:
         "external/mito-genotyping/exampleProcessing/02_merge_pileup_counts.sh {input.scPileup_dir} {params.samplename}"


rule scPileup_MT_matrix:
    """Create the position-by-cell coverage matrix"""
    input:
        all = "data/processed/{sample}/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.txt.gz",
        scPileup_dir = "data/processed/{sample}/scPileup_concat_{num_read}",
        barcode_p = "data/processed/{sample}/{sample}_barcode_data.p"
    output:
        sc_coverage_f = "data/processed/{sample}/{sample}_{num_read}/sc_coverage.csv"
    shell:
        "python src/plot_heatmap_coverage.py sc_mt {input.barcode_p} {input.scPileup_dir} {output.sc_coverage_f} {maxBP}"


rule plot_scPileup_MT_matrix:
    """Plot the posiitonal coverages."""
    input:
        sc_coverage_f = "data/processed/{sample}/{sample}_{num_read}/sc_coverage.csv"
    output:
        save_f_heat = "figures/{sample}/{sample}_{num_read}_MT_position.png",
        save_f_coverage = "figures/{sample}/{sample}_{num_read}_MT_position_coverage.png"
    shell:
        "python src/plot_heatmap_coverage.py plot {input.sc_coverage_f} {output.save_f_heat} {output.save_f_coverage}"


# rule filter_cells:
#"""Create a text file of cell barcodes to keep"""
#     input:
#     output:
#     shell:
#
# rule filter_positions:
#""" Create a text file of variants to keep
#     input:
#     output:
#     shell:
#
# rule generate_allele_frequencies:
#"""Create the AF-by-cell csv file"""
#     input:
#     output:
#     shell:


#rule plot_allele_frequencies:
#"""Plot the AF-by-cell heatmap"""
