from os.path import dirname
# rule all:
#     input:
#         expand("{sample}/mgatk/vireoIn/cellSNP.tag.AD.mtx", sample=config["samples"])
#

rule get_refAllele:
    #input: config["mt_ref_fa"],
    params: config['mgatk']["chrM_refAllele"]
    output: "{outdir}/chrM_refAllele.txt"
    shell: 'cp {params} {output}'


rule mgatk:
    """ Run both toSeurat and call variants in one script"""
    input:
        all = "{outdir}/{sample}.coverage.txt",
        refAllele = "{outdir}/chrM_refAllele.txt"
    output:
        vars_f = "{outdir}/mgatk/{sample}.variant.rds",
        vars_qc = report("{outdir}/mgatk/{sample}.variantQC.png", category="Variants")
    params:
        data_dir = lambda wildcards, input: dirname(input.all),
        sample = lambda wildcards: wildcards.sample,
        outdir = "mgatk"
        #nCell = lambda wildcards: wildcards.nCell
    shell:
        "./R_scripts/wrap_mgatk.R {params.data_dir} {params.outdir}/{params.sample} FALSE "


rule mgatk_to_vireoIn:
    input: "{outdir}/mgatk/{sample}.variant.rds"
    output: "{outdir}/mgatk/{sample}/vireoIn/cellSNP.tag.AD.mtx"
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        sample = lambda wildcards: wildcards.sample
    shell:
        "python src/mgatk_to_vireo.py {params.indir} {params.outdir} {params.sample}"

