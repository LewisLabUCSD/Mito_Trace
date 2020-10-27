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
#print(pd.read_table(config["samples"], dtype=str,sep=','))
samples = pd.read_table(config["samples"], dtype=str,sep=',').set_index(["sample"], drop=False)
#(samples)

#print('index',samples.index)
res = config["results"]
mq = config["mapq"]


#workdir: config["work_dir"]

#
# wildcard_constraints:
#     mapq="\d+",
#     cellr='True|False'


rule all:
    input:
        expand("{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_{num_read}/af_filt_{min_cells}_{min_reads}_{top_pos}_{top_cells}_{cell_mt_coverage}.p", results=res, sample=samples["sample"].values,
               mapq = mq, cellr_bc=config["use_cellr_barcode"],
               num_read=config["num_reads_filter"],
               min_cells=config["min_cells"],
               min_reads=config["min_reads"],
               top_pos = config["top_pos"],
               top_cells = config["top_cells"],
               cell_mt_coverage=config["cell_mt_coverage"]),


# def get_sample_bai(wildcards):
#     bai = os.path.join(samples.loc[wildcards.sample, "raw"],samples.loc[wildcards.sample, "bam"]) + ".bam.bai"
#     return bai


def get_sample_bam(wildcards):
    bam = samples.loc[wildcards.sample,"bam_f"]
    #print(bam)
    return bam # + ".bam"



rule simulate_single_clone:
    input:

    output:

# rule cp_bam:
#     """Move bam file to current location"""
#     input: get_sample_bam
#     output:
#         bam = "{results}/{sample}/00_bam/{sample}.bam",
#         #bai = "{results}/{sample}/00_bam/{sample}.bam.bai"
#     shell: "cp {input} {output}"
# #    run:
#         # shell("ln -s {input} {output.bam}"),
#         # shell("ln -s {input}.bai {output.bai}")
#

rule simulate_sensitivity:
    aggregate()