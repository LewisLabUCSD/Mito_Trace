from os.path import join, dirname, basename
from src.config import ROOT_DIR
from src.utils.parse_config import read_config_file, write_config_file
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

peaks=config["gatk_bed_regions"].keys()

mut_filt = params["som_vars"]["mutect"]["params"]

rule all:
  input:
    expand("{outdir}/regions_{peaks}/gatk_mutect/variants.vcf.gz",
           outdir=join(res, "somatic_variants"), peaks=peaks),
    #expand("{outdir}/all_variants.vcf", outdir=join(res, "somatic_variants")),
    #expand("{outdir}/preproc/merge_bam/pooled.sorted.bam",outdir=join(res, "somatic_variants")),
    expand("{outdir}/regions_{peaks}/merge_filt_bam/pooled.sorted.bam",
           outdir=join(res, "somatic_variants"), peaks=peaks),
    expand("{outdir}/regions_{peaks}/gatk_mutect/post/variants.filt_af.{af_low}_{af_high}_ad.{ad_low}.vcf.gz",
           outdir=join(res, "somatic_variants"), peaks=peaks,
           af_low=mut_filt["af_low"], af_high=mut_filt["af_high"],
           ad_low=mut_filt["ad_low"]),
    expand("{outdir}/regions_{peaks}/gatk_mutect/post/variants.varscan.filt.vcf.gz",
       outdir=join(res, "somatic_variants"), peaks=peaks)
    #expand("{outdir}/varCA_pt1/callers", outdir=join(res,"somatic_variants")),
    #expand("{outdir}/gatk/variants.vcf.gz",  outdir=join(res,"somatic_variants")),
    #expand("{outdir}/gatk_mutect/variants.vcf.gz", outdir=join(res,"somatic_variants")),


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
        temp("{outdir}/regions_{peaks}/preproc/bam/{s}.peaks.bam")
    params:
        bed_in = lambda wildcards: config["bed_regions"][wildcards.peaks],
        #bed_in = config["merged_regions_f"]
    shell: "bedtools intersect -wa -a {input.bam_f} -b {params.bed_in} > {output}"


rule merge_filt_bams:
    input:
        bam = expand("{{outdir}}/regions_{{peaks}}/preproc/bam/{s}.peaks.bam", s=samples.index),
        #bai = expand("{{outdir}}/preproc/{s}.peaks.bai", sample=samples.index)
    output:
        bam= temp("{outdir}/regions_{peaks}/merge_filt_bam/pooled.bam"),
        sort_bam = "{outdir}/regions_{peaks}/merge_filt_bam/pooled.sorted.bam"
    shell: "samtools merge {output.bam} {input.bam} && samtools sort {output.bam} > {output.sort_bam} && samtools index {output.sort_bam} "



rule genome_dict:
    input: ref_fa
    output: ref_fa.replace(".fa", "") + ".dict"
    shell: "gatk CreateSequenceDictionary -R {input}"



rule run_mutect_gatk:
    input:
        bam= "{outdir}/regions_{peaks}/merge_filt_bam/pooled.sorted.bam",
        ref_dict = ref_fa.replace(".fa", "") + ".dict"
        #regions =
    output:
        vcf="{outdir}/regions_{peaks}/gatk_mutect/variants.vcf.gz"
    params:
        genome = ref_fa
    shell: "gatk Mutect2 -R {params.genome} -I {input.bam} -O {output.vcf}"


rule filter_mutect:
    input:
        vcf="{outdir}/regions_{peaks}/gatk_mutect/variants.vcf.gz"
    output:
        vcf="{outdir}/regions_{peaks}/gatk_mutect/post/variants.filt_af.{af_low}_{af_high}_ad.{ad_low}.vcf.gz"
    params:
        genome = ref_fa
    shell:
        "gatk FilterMutectCalls --min-allele-fraction 0.05 -R {params.genome}  -V {input.vcf}  -O {output.vcf}"

rule unzip:
    input: "{outdir}/regions_{peaks}/gatk_mutect/variants.vcf.gz"
    output: "{outdir}/regions_{peaks}/gatk_mutect/variants.vcf"
    shell: "gunzip -c {input} > {output}"


rule filter_varscan_mutect:
    input:
        vcf="{outdir}/regions_{peaks}/gatk_mutect/variants.vcf"
    output:
        vcf="{outdir}/regions_{peaks}/gatk_mutect/post/variants.varscan.filt.vcf.gz"
    params:
        genome = ref_fa
    shell: "varscan somaticFilter {input.vcf} --output-file {output.vcf} --min-coverage 20 --min-avg-qual 25 --min-var-freq 0.05 --min-reads2 10 "


rule run_gatk:
    input:
        bam= "{outdir}/regions_{peaks}/merge_filt_bam/pooled.sorted.bam",
        ref_dict = ref_fa.replace(".fa", "") + ".dict"
        #regions =
    output:
        tmp="{outdir}/regions_{peaks}/gatk_haplotype/variants.vcf.gz",
        vcf="{outdir}/regions_{peaks}/gatk/variants.vcf.gz"
    params:
        genome = ref_fa
    shell:
         """
            gatk --java-options "-Xmx4g" HaplotypeCaller -R {params.genome} -I {input.bam}  -O {output.tmp}
            && \
            gatk --java-options "-Xmx4g" GenotypeGVCFs  --include-non-variant-sites  -R {params.genome} -V {output.tmp} -O {output.vcf}
        """
         # # -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation


def get_varca_cfg():
    if "varCA_path" in config:
        return join(config["varCA_path"], "configs", "config.yaml")
    else:
        return join(ROOT_DIR, "software/varCA/configs/config.yaml")

def get_varca_prep():
    if "varCA_path" in config:
        return join(config["varCA_path"], "configs", "prepare.yaml")
    else:
        return join(ROOT_DIR, "software/varCA/configs/prepare.yaml")

def get_snp_call():
    if "varCA" in params and "params" in params["varCA"]:
        return params["varCA"]["params"]["snp_callers"]
    return ["gatk-snp"]

def get_indel_call():
    if "varCA" in params and "params" in params["varCA"]:
        return params["varCA"]["params"]["indel_callers"]
    return ["gatk-indel"]


##snp_callers: [gatk-snp, varscan-snp, vardict-snp]
#indel_callers: [gatk-indel, varscan-indel, vardict-indel, pindel, illumina-strelka]
rule create_varca_config:
    """Create the input config for varCA to run"""
    input:
        get_varca_cfg(),  #join(ROOT_DIR, "software/varCA/configs/config.yaml"),
        "{outdir}/regions_{peaks}/merge_filt_bam/pooled.sorted.bam"
    output:
        cfg="{outdir}/regions_{peaks}/varCA/config.yaml",
        samples="{outdir}/regions_{peaks}/varCA/samples.tsv"
    params:
        snp_call = get_snp_call(),
        indel_call = get_indel_call()
    run:
        cfg = read_config_file(input[0])
        cfg["genome"] = ref_fa
        cfg["out"] = dirname(output['cfg'])
        cfg["sample_file"] = output['samples']
        write_config_file(output.cfg, cfg)
        samples = f"merged\t{output['samples']}"
        with open(output['samples'],'w') as f:
            f.write(samples)

rule run_varCA:
    """Run GATK with varCA. Setup bam files to run"""
    input:
        "{outdir}/regions_{peaks}/merge_filt_bam/pooled.sorted.bam",
        cfg="{outdir}/regions_{peaks}/varCA/config.yaml",
    output:
        "{outdir}/regions_{peaks}/varCA/all_variants.vcf"
    shell: "varCA {input}"


rule create_varca_prepare_config:
    """Create the input config for varCA to run"""
    input:
        get_varca_prep(),  #join(ROOT_DIR, "software/varCA/configs/config.yaml"),
        bam="{outdir}/regions_{peaks}/merge_filt_bam/pooled.sorted.bam"
    output:
        cfg="{outdir}/regions_{peaks}/varCA_pt1/prepare.yaml",
        samples="{outdir}/regions_{peaks}/varCA_pt1/samples.tsv"
    params:
        snp_call = get_snp_call(),
        indel_call = get_indel_call()
    run:
        cfg = read_config_file(input[0])
        cfg["genome"] = ref_fa
        cfg["out"] = dirname(output['cfg'])
        cfg["sample_file"] = output['samples']
        cfg["SAMP_NAMES"] = False
        cfg["caller_path"] = join(ROOT_DIR, "software/varCA")
        write_config_file(output.cfg, cfg)
        samples = f"merged\t{input['bam']}"
        with open(output['samples'],'w') as f:
            f.write(samples)


rule run_varca_prepare:
    input:
        cfg = "{outdir}/regions_{peaks}/varCA_pt1/prepare.yaml",
    output:
        directory("{outdir}/regions_{peaks}/varCA_pt1/callers")
    params:
        prep_rule = join(ROOT_DIR, "software", "varCA", "rules", "prepare.smk")
    shell: "snakemake -s {params.prep_rule} --use-conda -j 16 --configfile {input.cfg}"


