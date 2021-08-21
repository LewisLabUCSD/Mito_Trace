from os.path import dirname


rule all:
    input: expand("{sample}/mgatk/vireoIn/cellSNP.tag.AD.mtx",
                  results=config["results"], sample=config["samples"])

rule get_refAllele:
    #input: config["mt_ref_fa"],
    params: config["chrM_refAllele"]
    output: "{sample}/chrM_refAllele.txt"
    shell: 'cp {params} {output}'


rule mgatk:
    """ Run both toSeurat and call variants in one script"""
    input:
        all = "coverage.txt",
        refAllele = "chrM_refAllele.txt"
    output:
        vars_f = "{sample}/mgatk/{sample}.variant.rds",
        vars_qc = report("mgatk/{sample}.variantQC.png")
    params:
        data_dir=lambda wildcards, input: dirname(input.all),
        sample = lambda wildcards: wildcards.sample,
    shell:
        "./R_scripts/wrap_mgatk.R {params.data_dir} mgatk/{params.sample} FALSE"


rule mgatk_to_vireoIn:
    input: "{sample}/mgatk/{sample}.variant.rds"
    output: "{sample}/mgatk/vireoIn/cellSNP.tag.AD.mtx"
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        sample = lambda wildcards: wildcards.sample
    shell:
        "python src/mgatk_to_vireo.py {params.indir} {params.outdir} {params.sample}"

