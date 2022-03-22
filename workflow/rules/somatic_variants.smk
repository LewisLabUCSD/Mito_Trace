from os.path import join, dirname, basename
from src.config import ROOT_DIR
from src.utils.parse_config import read_config_file
import os
import numpy as np
from os.path import join, dirname
import pandas as pd
from snakemake.utils import min_version
from icecream import ic
min_version("6.0")
print('config', config)


########################################################################
# Setup parameters and outdir
########################################################################
res = join(config["outdir"], "pipeline", config["prefix"])

params = read_config_file(config["config"])
samples = pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)
anno_res = join(config["outdir"], "annotation", "data", config["prefix"])
ref_fa = params["genome_path"][config["genome"]] ["ref_fa"] #data/processed/genomes/mtMasked/GRCh38_MT_blacklist_A2_2020/fasta/genome.fa
#mt_ref_fa = params["genome_path"][config["genome"]]["mt_ref_fa"]
print('samples', samples.index)

rule all:
  input:
    #expand("{outdir}/all_variants.vcf", outdir=join(res, "somatic_variants")),
    expand("{outdir}/preproc/merge_bam/pooled.sorted.bam",outdir=join(res, "somatic_variants")),

rule move_bam:
    input:
        bam_files = expand("{mtscATAC_dir}/{s}/outs/possorted_bam.bam", mtscATAC_dir=config['mtscATAC_OUTDIR'], s=samples.index),
    output:
        bam_in = expand("{{outdir}}/preproc/bam/{s}.bam", s=samples.index),
        bai_in = expand("{{outdir}}/preproc/bam/{s}.bam.bai", s=samples.index)
    params:
        bam_indir = lambda wildcards, output: dirname(output.bam_in[0])
    run:
        for ind, curr_bam in enumerate(input.bam_files):
          cmd = f"ln -s {curr_bam} {output.bam_in[ind]}"
          print('cmd', cmd)
          os.system(cmd)
          cmd2 = f"cp {curr_bam+'.bai'} {output.bam_in[ind]+'.bai'}"
          print('cmd2', cmd2)
          os.system(cmd2)

#
# # peak files from the merged peaks in annotation
# rule needlestack:
#     input:
#         bam_in = expand("{{outdir}}/preproc/bam/{s}.bam", s=samples.index),
#         bai_in = expand("{{outdir}}/preproc/bam/{s}.bam.bai", s=samples.index),
#         bed_in = config["merged_peaks_needle_f"]  # "{outdir}/bed/peaks.bed"
#     output:
#         vcf = "{outdir}/all_variants.vcf"
#     params:
#         bam_indir = lambda wildcards, input: dirname(input.bam_in[0]),
#         mt_ref = ref_fa
#     shell:
#         "nextflow run iarcbioinfo/needlestack --bed {input.bed_in} --input_bams {params.bam_indir}/ --ref {params.mt_ref} --output_vcf {output.vcf}"


rule barcode_addnames:
    input:
        barcode_files = expand("{mtscATAC_dir}/{s}/outs/filtered_peak_bc_matrix/barcodes.tsv",
            mtscATAC_dir = config['mtscATAC_OUTDIR'], s=samples.index),

    output:
        barcode_files = expand("{{outdir}}/preproc/barcodes/{s}.barcodes.tsv", s=samples.index)
    #run:



rule filt_bams:
    input:
        bam_f = "{outdir}/preproc/bam/{s}.bam",
    output:
        "{outdir}/preproc/bam/{s}.peaks.bam"
    params:
        bed_in = config["merged_peaks_small_f"]
    shell: "bedtools intersect -wa -a {input.bam_f} -b {params.bed_in} > {output}"


rule merge_bams_v01:
    input:
        bam_in = expand("{{outdir}}/preproc/bam/{s}.peaks.bam", s=samples.index),
        barcode_files = expand("{mtscATAC_dir}/{s}/outs/filtered_peak_bc_matrix/barcodes.tsv",
                                mtscATAC_dir = [config['mtscATAC_OUTDIR']], s=samples.index),
        #bed_in = config["merged_peaks_small_f"]
        #barcode_files = expand("{{outdir}}/preproc/barcodes/{s}.barcodes.tsv", s=samples.index),
        #bam_files = expand("{mtscATAC_dir}/{s}/outs/possorted_bam.bam", mtscATAC_dir=config['mtscATAC_OUTDIR'], s=samples.index),
    output:
        #outdir = directory("{outdir}/preproc/merge_bam"),
        #bam = "{outdir}/preproc/merge_bam/pooled.sorted.bam",
        bam = expand("{{outdir}}/preproc/merge_bam/pooled.{s}.bam", s=samples.index)
    params:
        outdir = lambda wildcards, output: dirname(output.bam[0]),
        barcode_files = lambda wildcards, input: ",".join(input.barcode_files),
        bam_files = lambda wildcards, input: ",".join(input.bam_in),
        seed = 42,
        samples = ",".join(samples.index)
    shell: "python workflow/notebooks/merge_bam/merge_bams.py --samFiles {params.bam_files} --barcodeFiles {params.barcode_files} --outDir {params.outdir} --samples {params.samples} --randomSEED {params.seed}"

rule merge_bams_v02:
    input:
        bam = "{outdir}/preproc/merge_bam/pooled.{sample}.bam"
    output:
        bai = "{outdir}/preproc/merge_bam/pooled.{sample}.bai"
    shell: "samtools index {input.bam}"

rule merge_bams_v03:
    input:
        bam = expand("{{outdir}}/preproc/merge_bam/pooled.{sample}.bam", sample=samples.index),
        bai = expand("{{outdir}}/preproc/merge_bam/pooled.{sample}.bai", sample=samples.index)
    output:
        "{outdir}/preproc/merge_bam/pooled.sorted.bam"
    shell: "samtools merge {input.bam} -o {output} && samtools sort {output}"
