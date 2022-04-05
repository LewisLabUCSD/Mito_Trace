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
from src.config import ROOT_DIR
import subprocess as subp
########################################################################
# Setup parameters and outdir
########################################################################
res = join(config["outdir"], "pipeline", config["prefix"], "somatic_variants")


params = read_config_file(config["config"])
samples = pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)


rule all:
    input:
        # expand("{outdir}/gatk_mutect/post/{s}/scPileupVars",
        #     outdir=join(res), s=samples.index),
        expand("{outdir}/merge_filt_bam/_split_rg.txt", outdir=res),
        # expand("{outdir}/gatk_mutect/post/{s}/concat_scPileupVars/pileup.tsv",
        #         outdir=join(res), s=samples.index,
        #        ),
        expand("{outdir}/gatk_mutect/post/{s}/cells_vars/vcfpad_1/af.pileup.tsv",
                outdir=res, s=samples.index),
        expand("{outdir}/gatk_mutect/post/variants.annotate.vcf", outdir=res)


#############################################
# Convert needlestack results to cell-by-vars
############################################
def get_vcf(wildcards):
    return f"{wildcards.outdir}/gatk_mutect/variants.vcf.gz"
    # if wildcards.regions == "all":
    #     return f"{wildcards.outdir}/gatk_mutect/variants.vcf.gz"
    #
    # return f"{wildcards.outdir}/gatk_mutect/{wildcards.regions}.variants.vcf"

#"{outdir}/gatk_mutect/variants.vcf.gz"


# Need to install vep data after mamba
#vep_install -a cf -s homo_sapiens -y GRCh38 -c data/processed/genomes/homo_sapiens/vep_GRCh38 --CONVERT

rule bgzip_vcf:
    input:
        vcf = get_vcf
    output: "{outdir}/gatk_mutect/post/variants.vcf.bgz"
    shell: "zcat {input} | bgzip -c > {output} && tabix {output}"

rule annotate_vcf:
    input:
        vcf = "{outdir}/gatk_mutect/post/variants.vcf.bgz" # get_vcf
    output: "{outdir}/gatk_mutect/post/variants.annotate.vcf"
    params:
        gene_f = params["genome_path"][config["genome"]]["gene_sorted_file"]+".bgz",
        cache = join(ROOT_DIR, "data", "processed", "genomes", "homo_sapiens/vep_GRCh38/" )
    #shell: "bedtools annotate -i {input} -files {params.gene_f} > {output}"
    shell: "vep -i {input.vcf} --gff {params.gene_f} --cache --dir {params.cache} --assembly GRCh38 --offline --force_overwrite --no_progress --output_file {output}"


rule intersect_vcf:
    input:
        bam= "{outdir}/merge_filt_bam/pooled.sorted.bam",
        vcf =  "{outdir}/gatk_mutect/variants.vcf.gz" #get_vcf #"{outdir}/gatk_mutect/variants.all.vcf" #vcf = "{outdir}/gatk_mutect/{peaks}.variants.vcf"
    output:
        bam = ("{outdir}/gatk_mutect/post/merged.vcf.bam"),
        bai = ("{outdir}/gatk_mutect/post/merged.vcf.bam.bai"),
    resources:
        mem_mb=20000
    shell:
        "bedtools intersect -wa -a {input.bam} -b {input.vcf}  > {output.bam} && samtools index {output.bam}"


rule split_by_rg:
    input:
        bam= "{outdir}/gatk_mutect/post/merged.vcf.bam",
    output:
        tmpf = "{outdir}/merge_filt_bam/_split_rg.txt",
        #bams =  temp(expand("merged.vcf_{n}.bam", n=np.arange(len(samples))))
        bams =  temp(expand("{{outdir}}/merge_filt_bam/merged.vcf_{n}.bam", n=[str(x) for x in np.arange(len(samples))]))
        #bams = expand("{{outdir}}/merge_filt_bam/{s}", s=samples.index) #samples.index )
    params:
        outdir = lambda wildcards, output: dirname(output.tmpf)

    shell: "cd {params.outdir} && samtools split {input} && touch {output}"

rule index:
    input:
        bam = "{outdir}/merge_filt_bam/merged.vcf_{n}.bam"
    output:
        bai = "{outdir}/merge_filt_bam/merged.vcf_{n}.bam.bai"
    shell: "samtools index {input.bam}"


print('samples', [str(x) for x in np.arange(len(samples))])

if "rg_dict" not in config:
    config["rg_dict"] = "None"

rule get_cond_name:
    input:
        bams = expand("{{outdir}}/merge_filt_bam/merged.vcf_{n}.bam", n=[str(x) for x in np.arange(len(samples))]),
        bais = expand("{{outdir}}/merge_filt_bam/merged.vcf_{n}.bam.bai", n=[str(x) for x in np.arange(len(samples))]),
    output:
        bams = temp(expand("{{outdir}}/gatk_mutect/post/merged.vcf_{s}.bam", s=samples.index))
    params:
        samples = samples.index.values,
        outdir = lambda wildcards, output: dirname(output.bams[0]),   #.bams[0]),
        cond_d = config.get("rg_dict", "None")
    run:
        import pysam
        conversion = {}
        conversion_bai = {}
        print('in bams', input.bams)
        if params.cond_d is not "None":
            cond_d = {}
            for s in params.cond_d.split(";"):
                cond_d[s.split(",")[0]] = s.split(",")[1]
        else:
            cond_d = {s:s for s in samples.index}

        for in_f in input.bams:
            print(in_f)
            isInF = False
            samfile = pysam.AlignmentFile(in_f, "rb")
            for read in (samfile.fetch(until_eof=True)):
                # barcode itr for current read
                #if "CB" in read.tags():
                if isInF:
                    break
                curr_rg = read.get_tag('RG').split(":")[0]
                print('samples', params.samples)
                print('curr_rg', curr_rg)
                for s in params.samples:
                    if s in curr_rg:
                        print(f' {s} here')
                        if cond_d is not "None":
                            conversion[in_f] = join(params.outdir, f"merged.vcf_{cond_d[s]}.bam")
                            conversion_bai[in_f+".bai"] = join(params.outdir, f"merged.vcf_{cond_d[s]}.bam.bai")
                        else:
                            conversion[in_f] = join(params.outdir, f"merged.vcf_{s}.bam") #s
                            conversion_bai[in_f+".bai"] = join(params.outdir, f"merged.vcf_{s}.bam.bai")
                        isInF = True
                        break
            samfile.close()
        for c in conversion:
            cmd = f"mv {c} {conversion[c]}"
            os.system(cmd)
            cmd = f"mv {c} {conversion_bai[c+'.bai']}"
            os.system(cmd)



# rule filter_barcodes:
#     input:
#         bam = "{outdir}/gatk_mutect/post/{s}.bam",
#         bai = "{outdir}/gatk_mutect/post/{s}.bam.bai",
#     output:
#         bam = temp("{outdir}/gatk_mutect/post/cellFilt_{s}.bam"),
#         bai = temp("{outdir}/gatk_mutect/post/cellFilt_{s}.bam.bai"),

def get_barcode(wildcards):
    return samples.loc[wildcards.s, "barcode_f"]
    #barcode_files = expand("{mtscATAC_dir}/{s}/outs/filtered_peak_bc_matrix/barcodes.tsv",
            #mtscATAC_dir = config['mtscATAC_OUTDIR'], s=samples.index),
    #return allCellBarcodes[wildcards.s]

rule filter_barcodes_01:
    input:
        bam = "{outdir}/gatk_mutect/post/merged.vcf_{s}.bam",
        cells= get_barcode  #"{output}/data/{sample}/MT/cellr_True/{sample}_barcode_data.txt"
        #"{output}/data/{sample}/MT/{sample}.MT.bam",
    output:
        sam = temp("{outdir}/gatk_mutect/post/filtered.{s}.sam"),
    params:
        cmd = lambda wildcards, input: f"{{ samtools view -H {input.bam} & samtools view {input.bam} | LC_ALL=C grep -F -f {input.cells}; }}"
    resources:
        mem_mb=80000
    run: shell("{params.cmd} > {output.sam} 2> {output.sam}.log")
    #shell:  "{{ samtools view -H {input[0]} & samtools view {input[0]} | LC_ALL=C grep -F -f {input.cells}; }} > {output.sam} 2> {output.sam}.log"


rule filter_barcodes_02:
    input:
        sam = "{outdir}/gatk_mutect/post/filtered.{s}.sam"
    output:
        bam = temp("{outdir}/gatk_mutect/post/filtered.{s}.bam") #temp("{output}/data/{sample}/MT/filtered.bam"),
    shell: "samtools view -b {input.sam} > {output.bam}"


rule sortCB:
    input:
        bam= "{outdir}/gatk_mutect/post/filtered.{s}.bam",
    output: "{outdir}/gatk_mutect/post/{s}/{s}.CB.bam"  #temp("{output}/data/{sample}/MT/{sample}.MT.CB.bam")
    shell: "samtools sort -t CB {input.bam} > {output}"


rule split_by_cb:
    input:
        bam = "{outdir}/gatk_mutect/post/{s}/{s}.CB.bam" #bam = "{outdir}/gatk_mutect/post/{s}.bam",
    output:
        scBam = temp(directory("{outdir}/gatk_mutect/post/{s}/scBam")),
        #cb_bam = "{outdir}/gatk_mutect/post/{s}.CB.bam"
    params:
        script = join(ROOT_DIR, "src/mtpreproc/split_by_CB.py") #script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/needlestack_convert_vars_to_cells.ipynb")
    resources:
        mem_mb = 20000
    shell:
        "python {params.script} {input.bam} {output.scBam}"


rule cb_to_pileup:
    input:
        "{outdir}/gatk_mutect/post/{s}/scBam"
    output:
        temp(directory("{outdir}/gatk_mutect/post/{s}/scPileupVars")),
    params:
        script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/scPileup_counts.py"),
        nproc = lambda wildcards: config["nproc"] if "nproc" in config else 24
    threads: 16
    shell: "python {params.script} {input} {output} --nproc {params.nproc}"


rule scPileup_concat:
    """ Concat the pileups for each separate cell pileup file
    """
    input:
        scPileup_dir = "{outdir}/gatk_mutect/post/{s}/scPileupVars"
    output:
        "{outdir}/gatk_mutect/post/{s}/concat_scPileupVars/pileup.tsv"
    threads: 28
    shell: "find {input} -type f -name *_pileups.bq.tsv -exec cat {{}} \\; > {output}"


rule pileup_to_cell_vars:
    input:
        pileup = "{outdir}/gatk_mutect/post/{s}/concat_scPileupVars/pileup.tsv",
        vcf = get_vcf
        #"{outdir}/gatk_mutect/post.variants.vcf"
    output:
        note = "{outdir}/gatk_mutect/post/{s}/cells_vars/vcfpad_1/out.ipynb",
        ref="{outdir}/gatk_mutect/post/{s}/cells_vars/vcfpad_1/af.ref.pileup.tsv",
        af="{outdir}/gatk_mutect/post/{s}/cells_vars/vcfpad_1/af.pileup.tsv"
    params:
        outdir = lambda wildcards, output: dirname(output.note),
        script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/pileup_to_af_dp.ipynb"),
        vcf_pad = 1
    resources:
        mem_mb=80000
    shell: "papermill -p vcf {input.vcf} -p scPileupVars {input.pileup} -p vcf_pad {params.vcf_pad} -p outdir {params.outdir} {params.script} {output.note}"



#cells_vars/vcfpad_1
