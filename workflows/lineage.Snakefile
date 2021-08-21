from src.config import ROOT_DIR
import pandas as pd


samples = pd.read_table(config["samples"], dtype=str,sep=',').set_index(["sample_name"], drop=False)

rule all:
    input:
        # 1. Lineage labels for cells, which can be input to loupe browser
        expand("clones/{clone_type}_{clone_type_params}/donor{d}/cells_BC.csv",
               clone_type=config["method"],
               clone_type_params=[config[x]['params'] for x in config["method"]],
               d=range(config["N_DONORS"]))
         ,

        # 2. Clones: types of clones
        expand("clones/{clone_type}_{clone_type_params}/variants.ipynb",
                clone_type=config["method"],
                clone_type_params=[config[x]['params'] for x in config["method"]],
                d=range(config["N_DONORS"])
               ),

        # 3. Enrichment: enrichment of each clone across donor
        expand("clones/{clone_type}_{clone_type_params}/enrichment/volcano_Fisher_foldNorm.png",
            clone_type=config["method"],
            clone_type_params=[config[x]['params'] for x in config["method"]])



########################################################################
## Workflow A: After multiplexing, get the union of the vcfs, and also
# run the variants for each donor
##
#rule clones_pipe:
    #"{results}/clones_pipe/{cloneMethod}/cells_meta.tsv"
from os.path import join, dirname
import os


def get_multiplex_input(wildcards):
    return config["files"]["multiplex"]


########################################################################
## Workflow A: Directly from multiplex output
########################################################################
rule vireo:
    input: config["files"]["multiplex"]
    output:
        report(expand("clones/variants_simple/{{clone_type}}_{clone_type_params}/donor{d}OUT.labels.png",
                      d=range(config["N_DONORS"]), n_clone=config['clone']['vireo']['params']['nclonelist'],
                      clone_type_params=config["simple"]["params"])),

        report(expand("clones/variants_simple/{{clone_type}}_{clone_type_params}/donor{d}OUT.variants.labels.png",
               d=range(config["N_DONORS"]), clone_type_params=config["simple"]["params"])),
    params:
        notebook=join("src", "vireo", "2_MT_Lineage_Construct.ipynb"),
        output_notebook = lambda wildcards, output: join(dirname(dirname(output[0])), 'clones.ipynb'),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["N_DONORS"],
        n_clone= ",".join([str(x) for x in config['clone']['vireo']['params']['nclonelist']])
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p n_clone_list {params.n_clone} {params.notebook} {params.output_notebook}"



########################################################################
## Intermediate step for other workflows.
## Extract pileups for each donor from original filter matrices and run mgatk
########################################################################
def get_coverage(wildcards):
    return f"{config['cov_indir']['sample']}/{wildcards.sample}.coverage.txt"

def get_cells_meta(wildcards):
    return f"{config['cov_indir']['sample']}/{wildcards.sample}.coverage.txt"

rule scPileup_filter_mgatk:
    input:
        cov = config['files']['filters']['coverage'].values(),
        mult = join("config['files']['multiplex']",".ipynb") #"clones/{clone_type}_{clone_type_params}/multiplex.ipynb"
    output:
        cov = expand("clones/donors_mgatk_in/donor{d}/d{d}.coverage.txt",
                     d=range(config["N_DONORS"]))
    params:
        #outdir = lambda wildcards, output: dirname(output.cov),
        cells_meta = lambda wildcards, input: join(dirname(input.mult), "cells_meta.tsv"),
        sample = samples['sample_name'].values
    script: join(ROOT_DIR, "src/donor_filter_mgatk.py")


rule donor_get_refAllele:
    #input: config["mt_ref_fa"],
    params: config["chrM_refAllele"]
    output: "clones/donors_mgatk_in/donor{d}/chrM_refAllele.txt"
    shell: 'cp {params} {output}'


rule donor_mgatk:
    input:
        "clones/donors_mgatk_in/donor{d}/d{d}.coverage.txt",
        rules.donor_get_refAllele.output
    output: "clones/donors_mgatk_in/donor{d}/donor_mgatk/d{d}.variant.rds"
    params:
        data_dir =lambda wildcards, input: dirname(input[0]),
        donor = lambda wildcards: f"d{wildcards.d}",
    shell:
        "./R_scripts/wrap_mgatk.R {params.data_dir} donor_mgatk/{params.donor} FALSE"

rule donor_mgatk_to_vireoIn:
    input: "clones/donors_mgatk_in/donor{d}/donor_mgatk/d{d}.variant.rds"
    output: "clones/variants_mgatkrerun/{clone_type}_{clone_type_params}/donor{d}/cellSNP.tag.AD.mtx",
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        donor = lambda wildcards: f"d{wildcards.d}"
    shell:
        "python src/mgatk_to_vireo.py {params.indir} {params.outdir} {params.donor}"


## Workflow B: After multiplexing, separate by donor, grab variants from filters
## that overlap with both conditions, and then run mgatk to call variants again.
##
# Extract pileups for each donor from original filter matrices
########################################################################
rule donor_copy_cells_meta:
    input: config['files']['multiplex']
    output:
        "clones/variants_mgatkrerun/{clone_type}_{clone_type_params}/donor{d}/cells_meta.tsv"
    params:
        cells_meta = lambda wildcards, input: join(dirname(input[0]), "cells_meta.tsv")
    shell: "cp {params.cells_meta} {output}"


rule vireo_donor_mgatk:
    input:
        "clones/variants_mgatkrerun/{clone_type}_{clone_type_params}/donor{d}/cellSNP.tag.AD.mtx",
        "clones/variants_mgatkrerun/{clone_type}_{clone_type_params}/donor{d}/cells_meta.tsv"
    output: "clones/variants_mgatkrerun/{clone_type}_{clone_type_params}/results/clones.ipynb"
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        donor=lambda wildcards: wildcards.d, #config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "2b_MT_Lineage_Construct_mgatkDonors.ipynb"),
        workdir = os.getcwd(),
        n_clone=",".join([str(x) for x in config['clone']['vireo']['params']['nclonelist']])
    shell: "papermill --cwd {params.workdir} -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p donor {params.donor} -p n_clone_list {params.n_clone} {params.notebook} {output}"



########################################################################
## Workflow C: Using the same MGATK variants as above but running KNN instead
########################################################################

########################################################################

########################################################################
## For each cluster workflow, prepare the input into downstream analysis
########################################################################
#rule clones_donor_mgatk_process:
def check_cluster_method(wildcards):
    if wildcards.clone_type == "simple":
        return f"clones/variants_simple/{wildcards.clone_type}_{wildcards.clone_type_params}/results/clones.ipynb"
    elif wildcards.clone_type == "mgatkrerun_vireo":
        return f"clones/variants_mgatkrerun/{wildcards.clone_type}_{wildcards.clone_type_params}/results/clones.ipynb"


rule clone_method:
    input:
        check_cluster_method
    output:
        "clones/{clone_type}_{clone_type_params}/results/clones.ipynb"
    shell:
        "cat {input} > {output}"


########################################################################
## Break up variants by clones and get the types of variants
rule clones_type_variants:
    input: "clones/{clone_type}_{clone_type_params}/results/clones.ipynb"
    output: "clones/{clone_type}_{clone_type_params}/variants.ipynb"
    params:
        notebook=join("src", "vireo", join("6_MT_Clones_variantTypes.ipynb")),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS = config["N_DONORS"],
        sample_names = ','.join(config["sample_dict"].keys()), # make it as a list
        n_clones = lambda wildcards: wildcards.n_clones, #config['multiplex']["n_clone_list"],#lambda wildcards: wildcards.n_clones,
        var_thresh=0.001,
        vars_to_plot=10
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR}  -p n_clones {params.n_clones} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} -p var_thresh {params.var_thresh} -p vars_to_plot {params.vars_to_plot} {params.notebook} {output} && jupyter nbconvert --to pdf {output}"


########################################################################
rule enrichment:
    input:
        expand("clones/{{clone_type}}/clones_{{clone_type_params}}/donor{d}OUT.labels.png",
                        d=range(config["N_DONORS"]), n_clone=config['clone']['vireo']['params']['nclonelist'])
    output:
        report(expand("clones/{{clone_type}}/clones_{{clone_type_params}}/enrichment/volcano_Fisher_foldNorm.png",
                        n_clone=config['clone']['vireo']['params']['nclonelist'])),
    params:
        clones_indir = lambda wildcards, input: dirname(dirname(dirname(input[0]))),#lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        n_clones = config['simple']["params"]["n_clone_list"],
        script = join("src", "lineage_enrichment.py"),
        samples=",".join(config["multiplex"]['samples'])
    shell: "python {params.script} {params.clones_indir} {params.OUTDIR} {params.n_clones} {params.samples}"

