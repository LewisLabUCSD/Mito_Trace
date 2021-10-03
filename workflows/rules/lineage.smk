from src.config import ROOT_DIR
import pandas as pd
from os.path import join, dirname
import os
import numpy as np

wildcard_constraints:
    variants = "simple|mgatkdonor",
    d = "%d"


samples = config["samples"] #pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)

clones_cfg = config["clones"]

nclonelist = clones_cfg['vireo']['params']['nclonelist']

def get_input_multiplex(wildcards):
    if "files" in config:
        return config["files"]["multiplex"]["notebook"]
    else:
        return f"{wildcards.outdir}/multiplex/multiplex.ipynb"


def get_input_coverage(wildcards):
    if "files" in config:
        return config["files"]["mtpreproc"]["coverage"][wildcards.sample]
    else:
        return f"{wildcards.outdir}/{wildcards.sample}.coverage.txt"
    return




########################################################################
## Variant Intermediate step for other workflows.
## Extract pileups for each donor from original filter matrices and run mgatk
########################################################################
def get_coverage(wildcards):
    return f"{config['cov_indir']['sample']}/{wildcards.sample}.coverage.txt"

def get_cells_meta(wildcards):
    return f"{config['cov_indir']['sample']}/{wildcards.sample}.coverage.txt"



rule scPileup_filter_mgatk:
    input:
        cov =  expand("{{output}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/{sample}.coverage.txt",
                      sample=samples.index),
        mult = "{outdir}/mgatk/vireoIn/multiplex/multiplex.ipynb",
    output:
        cov = expand("{{outdir}}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/d{d}.coverage.txt",
                     d=range(config["N_DONORS"]))
    params:
        #outdir = lambda wildcards, output: dirname(output.cov),
        cells_meta = lambda wildcards, input: join(dirname(input.mult), "cells_meta.tsv"),
        sample = samples['sample_name'].values
    script: join(ROOT_DIR, "src/donor_filter_mgatk.py")


# rule donor_get_refAllele:
#     #input: config["mt_ref_fa"],
#     params: params["mgatk"]["chrM_refAllele"]
#     output: "{outdir}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/chrM_refAllele.txt"
#     shell: 'cp {params} {output}'
#

rule donor_mgatk:
    input:
        "{outdir}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/d{d}.coverage.txt",
        "{outdir}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/chrM_refAllele.txt"
    output: "{outdir}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/donor_mgatk/d{d}.variant.rds"
    params:
        data_dir =lambda wildcards, input: dirname(input[0]),
        donor = lambda wildcards: f"d{wildcards.d}",
    shell:
        "./R_scripts/wrap_mgatk.R {params.data_dir} donor_mgatk/{params.donor} FALSE"


########################################################################
## Workflow A: Directly from multiplex output
########################################################################
rule vireo:
    input:
        "{outdir}/mgatk/vireoIn/multiplex/multiplex.ipynb" #get_input_multiplex
    output:
        report("{outdir}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{nclones}/donor{d}.labels.png"),
        report("{outdir}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{nclones}/donor{d}.variants.labels.png")
        # report(expand("{{outdir}}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{{nclones}}/donor{d}.labels.png",d=np.arange(config["N_DONORS"]))),
        # report(expand("{{outdir}}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{{nclones}}/donor{d}.variants.labels.png",d=np.arange(config["N_DONORS"]))),
        #"{outdir}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{nclones}/donor{d}_clones.ipynb", #"{outdir}/mgatk/vireoIn/clones/variants_simple/vireo/clones.ipynb",
    params:
        notebook=join("src", "vireo", "2_MT_Lineage_Construct.ipynb"),
        #output_notebook = lambda wildcards, output: join(dirname(output[0]), 'clones.ipynb'),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["N_DONORS"],
        nclones= lambda wildcards: wildcards.nclones #",".join([str(x) for x in nclonelist])
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p nclones {params.nclones} {params.notebook} {params.OUTDIR}/output.ipynb"


########################################################################
## Workflow B: After multiplexing, separate by donor, grab variants from filters
## that overlap with both conditions, and then run mgatk to call variants again.
##
# Extract pileups for each donor from original filter matrices
########################################################################
rule donor_mgatk_to_vireoIn:
    input: "{outdir}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/donor_mgatk/d{d}.variant.rds"
    output: "{outdir}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/donor{d}/cellSNP.tag.AD.mtx",
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        donor = lambda wildcards: f"d{wildcards.d}"
    shell:
        "python src/mgatk_to_vireo.py {params.indir} {params.outdir} {params.donor}"


# def get_input_multiplex_cells(wildcards):
#     if "files" in config:
#         return config["files"]["multiplex"]["cells_meta"]
#     else:
#         return f"{wildcards.outdir}/multiplex/cells_meta.tsv"

rule donor_copy_cells_meta:
    input: "{outdir}/mgatk/vireoIn/multiplex/multiplex.ipynb"
    output:
        "{outdir}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/donor{d}/cells_meta.tsv"
    params:
        cells_meta = lambda wildcards, input: join(dirname(input[0]), "cells_meta.tsv")
    shell: "cp {params.cells_meta} {output}"


rule vireo_donor_mgatk:
    input:
        mtx=expand("{{outdir}}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/donor{d}/cellSNP.tag.AD.mtx",
                   d=np.arange(config["N_DONORS"])),
        cells=expand("{{outdir}}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/donor{d}/cells_meta.tsv", d=np.arange(config["N_DONORS"])),
    #output: "clones/variants_mgatkdonor/vireo/clones.ipynb"
    output:
        #"{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/clones.ipynb",
        report("{outdir}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}.labels.png"),
        report("{outdir}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}.variants.labels.png"),
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        donor=lambda wildcards: wildcards.d, #config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "2b_MT_Lineage_Construct_mgatkDonors.ipynb"),
        workdir = os.getcwd(),
        nclones= lambda wildcards: wildcards.nclones
    shell: "papermill --cwd {params.workdir} -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p donor {params.donor} -p nclones {params.nclones} {params.notebook} {params.OUTDIR}/output.ipynb"

# Need extra rule to make it like it's finished

########################################################################
## Break up variants by clones and get the types of variants
rule clones_type_variants:
    input:
        expand("{{outdir}}/mgatk/vireoIn/clones/variants_{{variants}}/vireo/nclones{{nclones}}/donor{d}.labels.png", d=np.arange(config["N_DONORS"]))
    output: "{outdir}/mgatk/vireoIn/clones/variants_{variants}/{method}/nclones{nclones}/variants.ipynb"
    params:
        notebook=join("src", "vireo", join("6_MT_Clones_variantTypes.ipynb")),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS = config["N_DONORS"],
        sample_names = ','.join(samples.index), # make it as a list
        nclones = lambda  wildcards: wildcards.nclones, #lambda wildcards: wildcards.nclones, #config['multiplex']["n_clone_list"],#lambda wildcards: wildcards.nclones,
        var_thresh=0.001,
        vars_to_plot=10
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR}  -p nclones {params.nclones} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} -p var_thresh {params.var_thresh} -p vars_to_plot {params.vars_to_plot} {params.notebook} {output} && jupyter nbconvert --to pdf {output}"


########################################################################
rule enrichment_vireo:
    input:
        expand("{{outdir}}/mgatk/vireoIn/clones/variants_{{variants}}/vireo/nclones{{nclones}}/donor{d}.labels.png", d=np.arange(config["N_DONORS"]))
    output:
        report("{outdir}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/enrichment/volcano_Fisher_foldNorm.png"),
    params:
        clones_indir = lambda wildcards, input: dirname(dirname(dirname(input[0]))),#lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        nclones = lambda wildcards: wildcards.nclones, #clones_cfg['vireo']["params"]["nclonelist"],
        script = join("src", "lineage", "lineage_enrichment.py"),
        samples=",".join(samples.index)
    shell: "python {params.script} {params.clones_indir} {params.OUTDIR} {params.nclones} {params.samples}"


# Bring enrichment results into same space
rule vireo_process:
    input:
        expand("{{outdir}}/mgatk/vireoIn/clones/variants_{{variants}}/vireo/nclones{nclones}/{f_out}",
                    zip,
                    nclones=nclonelist, f_out=["enrichment/volcano_Fisher_foldNorm.png", "variants.ipynb"]),
        #expand("{{outdir}}/mgatk/vireoIn/clones/variants_{{variants}}/vireo/nclones{nclones}/
    output:
        temp("{outdir}/mgatk/vireoIn/clones/variants_{variants}/vireo/temp/.tmp"),
    shell: "touch {output}"


# Sort of like an all rule to bring downstream results together.
rule complete_lineage:
    input:
        "{outdir}/mgatk/vireoIn/clones/variants_{variants}/{method}/temp/.tmp",
    output: "{outdir}/mgatk/vireoIn/clones/variants_{variants}/{method}/.completed"
    shell: "touch {output}"
