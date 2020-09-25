wildcard_constraints:
    num_read="\d+"

import os
import pandas as pd
from snakemake.utils import validate

#configfile: "parameters/Zhu_Single_Cell_10X_genomics_Human2_002.yaml"

#raw = config["raw"]

#bam = config["bam"]
num_reads_filter = config["num_reads_filter"]
maxBP = config["maxBP"]
ref_fa = config["ref_fa"]
print(pd.read_table(config["samples"], dtype=str,sep=','))
samples = pd.read_table(config["samples"], dtype=str,sep=',').set_index(["sample"], drop=False)
#(samples)

#print('index',samples.index)


#workdir: config["work_dir"]


wildcard_constraints:
    mapq="\d+",
    cellr='True|False'

rule all:
    input:
        #expand("{results}/{sample}/mapq_{mapq}/{sample}_scPileup_{num_read}",results=config["results"],sample=samples["sample"].values, num_read=config["num_reads_filter"], mapq=config["mapq"]),
        expand("{results}/{sample}/mapq_{mapq}/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.txt.gz",results=config["results"],sample=samples["sample"].values, num_read=config["num_reads_filter"], mapq=config["mapq"]),
        expand("{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_CB_coverage_hist_minReads{num_read}.png",results=config['results'],sample=samples["sample"].values, mapq=config["mapq"], cellr_bc=config["use_cellr_barcode"], num_read=config["num_reads_filter"]),
        expand("{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_{num_read}_MT_position.png",results=config['results'], sample=samples["sample"].values, num_read=config["num_reads_filter"], mapq=config["mapq"], cellr_bc=config["use_cellr_barcode"]),
        expand("{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_{num_read}_MT_position_coverage.png", results=config['results'],sample=samples["sample"].values, num_read=config["num_reads_filter"], mapq=config["mapq"], cellr_bc=config["use_cellr_barcode"])



# def get_sample_bai(wildcards):
#     bai = os.path.join(samples.loc[wildcards.sample, "raw"],samples.loc[wildcards.sample, "bam"]) + ".bam.bai"
#     return bai


def get_sample_bam(wildcards):
    bam = samples.loc[wildcards.sample,"bam_f"]
    #print(bam)
    return bam # + ".bam"


def get_sample_barcodes(wildcards):
    return os.path.join(os.path.dirname(samples.loc[wildcards.sample,"bam_f"]), "filtered_feature_bc_matrix/barcodes.tsv.gz")


def results_dir():
    return config["results"]

# rule orig_index:
#     input:  get_raw_bam #"{raw_f}"  #lambda wildcards: f"{config['bam'][wildcards.sample]}"
#     output: "{bam_f}.bam.bai" #"{raw}/{bam}.bai"
#     shell: "samtools index {input}"


rule cp_bam:
    """Move bam file to current location"""
    input: get_sample_bam
    output:
        bam = "{results}/{sample}/00_bam/{sample}.bam",
        #bai = "{results}/{sample}/00_bam/{sample}.bam.bai"
    shell: "cp {input} {output}"
#    run:
        # shell("ln -s {input} {output.bam}"),
        # shell("ln -s {input}.bai {output.bai}")


rule index_bam:
    """Index the bam file"""
    input: "{results}/{sample}/00_bam/{sample}.bam"
    output: "{results}/{sample}/00_bam/{sample}.bam.bai"
    shell: "samtools index {input}"


rule filter_bam:
    """Filter reads with low MAPQ and other potential filters """
    input:
        in_bam = "{results}/{sample}/00_bam/{sample}.bam",
        in_bai = "{results}/{sample}/00_bam/{sample}.bam.bai",
    output:
        mq_bam = "{results}/{sample}/01_bam_filter/mapq_{mapq}/{sample}.bam",
        mq_bai = "{results}/{sample}/01_bam_filter/mapq_{mapq}/{sample}.bam.bai"
    run:
            shell("samtools view {input.in_bam} -q {wildcards.mapq} -h -b > {output.mq_bam}")
            shell("samtools index {output.mq_bam}")


rule MT_map:
    """Extract the MT genome"""
    input:
        bam = "{results}/{sample}/01_bam_filter/mapq_{mapq}/{sample}.bam",
        bai = "{results}/{sample}/01_bam_filter/mapq_{mapq}/{sample}.bam.bai"
    output:
        mt_bam="{results}/{sample}/mapq_{mapq}/{sample}.MT.bam",
        mt_bai="{results}/{sample}/mapq_{mapq}/{sample}.MT.bam.bai"
    run:
        #shell("samtools {input.bam}")
        shell("samtools view -b {input.bam} MT > {output.mt_bam}")
        shell("samtools index {output.mt_bam}")


rule barcode_data:
    """Loop through the bam file and extract the barcode information."""
    input:
        mt_bam="{results}/{sample}/mapq_{mapq}/{sample}.MT.bam",
        mt_bai="{results}/{sample}/mapq_{mapq}/{sample}.MT.bam.bai"
    output: "{results}/{sample}/mapq_{mapq}/{sample}_barcode_data.p",
    shell:
         "python src/bam_barcodes_function.py {input.mt_bam} {output}"


rule barcode_filter:
    input:
        barcode_p = "{results}/{sample}/mapq_{mapq}/{sample}_barcode_data.p",
        cellr_f = get_sample_barcodes
    output: "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_barcode_data.p"
    params: cellr_bc = "{cellr_bc}"
    shell: "python src/filter_barcodes.py {input} {params} {output}"


## TODO
# rule compare_CB_with_cellranger_CB_list:
#     input:
#         cb_list = "{results}/{sample}/mapq_{mapq}/{sample}_barcode_data.p",
#         cellranger_list = "{results}/{sample}/mapq_{mapq}/{sample}.barcodes.txt"
#     output:
#         "{results}/{sample}/mapq_{mapq}/barcodes_compare/barcodes_detected.png"
#     shell:
#

# rule plot_CB_coverage:
#     """Plot the MT coverage of the single cells"""
#     input: "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_barcode_data.p"
#     output: "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_CB_coverage_hist.png"
#     shell:
#         "python src/plot_CB_coverage.py {input} {output}"


rule scBam:
    """Extract each single-cell and put into respective bam file"""
    input:
        mt_bam="{results}/{sample}/mapq_{mapq}/{sample}.MT.bam",
        mt_bai="{results}/{sample}/mapq_{mapq}/{sample}.MT.bam.bai"
    output: directory("/data/isshamie/mito_lineage/{results}/{sample}/mapq_{mapq}/{sample}_scBam")
    #directory("/data/isshamie/miito_lineage/mttrace/{sample}/mapq_{mapq}/{sample}_scBam")
    shell:
        "python src/split_by_CB.py {input.mt_bam} {output}"


rule scPileup:
    """Run the first part of the MT-genotype function by getting the read pileups for each bam file for each nucleotide and overall coverage"""
    input:
        scBam = "/data/isshamie/mito_lineage/{results}/{sample}/mapq_{mapq}/{sample}_scBam",
        #scBam = "{results}/{sample}/mapq_{mapq}/{sample}_scBam",
        barcodes = "{results}/{sample}/mapq_{mapq}/{sample}_barcode_data.p",
    output:
          directory("{results}/{sample}/mapq_{mapq}/{sample}_scPileup_{num_read}")
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
        scPileup_dir = "{results}/{sample}/mapq_{mapq}/{sample}_scPileup_{num_read}",
        # Added to save directly to /data/isshamie
        scBam = "/data/isshamie/mito_lineage/{results}/{sample}/mapq_{mapq}/{sample}_scBam"
    output:
        all = "{results}/{sample}/mapq_{mapq}/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.txt.gz"
    params:
        samplename = lambda wildcards, output: output.all.split("_all.coverage.txt.gz")[0]
          #"data/processed/{sample}/scPileup_concat/{sample}_{num_read}"
    shell:
         "software/mito-genotyping/exampleProcessing/02_merge_pileup_counts.sh {input.scPileup_dir} {params.samplename}"



rule plot_CB_coverage:
    """Plot the MT coverage of the single cells"""
    input:
         all = "{results}/{sample}/mapq_{mapq}/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.txt.gz",
         barcodes = "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_barcode_data.p"
    output: "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_CB_coverage_hist_minReads{num_read}.png"
    shell:
        "python src/plot_CB_coverage.py {input} {output}"


rule scPileup_MT_matrix:
    """Create the position-by-cell coverage matrix"""
    input:
        all = "{results}/{sample}/mapq_{mapq}/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.txt.gz",
        scPileup_dir = "{results}/{sample}/mapq_{mapq}/{sample}_scPileup_{num_read}",
        barcode_p = "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_barcode_data.p"
    output:
        sc_coverage_f = "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_{num_read}/sc_coverage.csv"
    shell:
        "python src/plot_heatmap_coverage.py sc_mt {input.barcode_p} {input.scPileup_dir} {output.sc_coverage_f} {maxBP}"


rule plot_scPileup_MT_matrix:
    """Plot the posiitonal coverages."""
    input:
        sc_coverage_f = "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_{num_read}/sc_coverage.csv"
    output:
        save_f_heat = "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_{num_read}_MT_position.png",
        save_f_coverage = "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_{num_read}_MT_position_coverage.png"
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
