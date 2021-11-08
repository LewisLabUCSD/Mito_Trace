from src.config import ROOT_DIR
import pandas as pd
from os.path import join, dirname
import os
import numpy as np
from icecream import ic


wildcard_constraints:
    variants = "simple|mgatkdonor",
    d = "%d"


samples = config["samples"] #pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)

clones_cfg = config["params"]["clones"]

nclonelist = clones_cfg['vireo']['params']['nclonelist']

# def get_input_multiplex(wildcards):
#     if "files" in config:
#         return config["files"]["multiplex"]["notebook"]
#     else:
#         return f"{wildcards.outdir}/multiplex/multiplex.ipynb"
#
#
# def get_input_coverage(wildcards):
#     if "files" in config:
#         return config["files"]["mtpreproc"]["coverage"][wildcards.sample]
#     else:
#         return f"{wildcards.outdir}/{wildcards.sample}.coverage.txt"
#     return
#
#


########################################################################
## Variant Intermediate step for other workflows.
## Extract pileups for each donor from original filter matrices and run mgatk
########################################################################
# def get_coverage(wildcards):
#     return f"{config['cov_indir']['sample']}/{wildcards.sample}.coverage.txt"

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





########################################################################
## Vireo A: Directly from multiplex output
########################################################################
rule vireo:
    input:
        "{outdir}/multiplex/multiplex.ipynb" #get_input_multiplex
    output:
        expand("{{outdir}}/variants_simple/vireo/nclones{{nclones}}/donor{d}.labels.png", d=np.arange(config["N_DONORS"])),#, category="lineage"),
        expand("{{outdir}}/variants_simple/vireo/nclones{{nclones}}/donor{d}.variants.labels.png", d=np.arange(config["N_DONORS"])),#, category="lineage"),
        "{outdir}/clones/variants_simple/vireo/nclones{nclones}/cells_meta.tsv"
    params:
        notebook=join("src", "vireo", "2_MT_Lineage_Construct.ipynb"),
        #output_notebook = lambda wildcards, output: join(dirname(output[0]), 'clones.ipynb'),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["N_DONORS"],
        nclones= lambda wildcards: wildcards.nclones #",".join([str(x) for x in nclonelist])
    threads: 8
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p nclones {params.nclones} {params.notebook} {params.OUTDIR}/output.ipynb"


########################################################################
## Vireo B: After multiplexing, separate by donor, grab variants from filters
## that overlap with both conditions, and then run mgatk to call variants again.
##
# Extract pileups for each donor from original filter matrices
########################################################################
rule vireo_mgatkdonor:
    input:
        mtx= "{outdir}/clones/variants_mgatkdonor/vireo/donor{d}/mgatk_donor/cellSNP.tag.AD.mtx",
        cells="{outdir}/clones/variants_mgatkdonor/vireo/donor{d}/mgatk_donor/cells_meta.tsv",
    #output: "clones/variants_mgatkdonor/vireo/clones.ipynb"
    output:
        #"{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/clones.ipynb",
        "{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}.labels.png",#, category="lineage", subcategory="donorMGATK"),
        "{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}.variants.labels.png",# category="lineage", subcategory="donorMGATK"),
        "{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}_cells_meta.tsv",
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        #N_DONORS =config["N_DONORS"], #config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "2b_MT_Lineage_Construct_mgatkDonors.ipynb"),
        d = lambda wildcards: wildcards.d,
        #workdir = os.getcwd(),
        nclones= lambda wildcards: wildcards.nclones
    threads: 8
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p donor {params.d} -p nclones {params.nclones} {params.notebook} {params.OUTDIR}/output.ipynb"


rule vireo_mgatkdonor_concat:
    input:
        expand("{{outdir}}/variants_mgatkdonor/vireo/nclones{{nclones}}/donor{d}_cells_meta.tsv",
               d = np.arange(config["N_DONORS"])),
    output:
         "{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/cells_meta.tsv",
    run:
        all = []
        for i in input:
            all.append(pd.read_csv(i, sep='\t'))
        all = pd.concat(all, ignore_index=True).sort_values(["donor", "lineage", "donor_index", "lineage_index"])
        if 'level_0' in all.columns:
                all = all.drop('level_0', axis=1)
        ic('all')
        ic(all.head())
        all.to_csv(output[0], sep='\t', index=False)

## Enrichment
rule vireo_enrichment:
    version: '1.0' # no +1 norm error
    input:
        "{outdir}/clones/variants_{variants}/vireo/nclones{nclones}/cells_meta.tsv"
    output:
        report("{outdir}/clones/variants_{variants}/vireo/nclones{nclones}/enrichment/volcano_Fisher_foldNorm.png", category="enrichment"),
    params:
        #clones_indir = lambda wildcards, input: dirname(input[0]),#lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        script = join("src", "lineage", "lineage_enrichment.py"),
        samples=",".join(samples.index)
    shell: "python {params.script} {input} {params.OUTDIR} {params.samples}"


# Bring enrichment results into same space
rule vireo_process:
    input:
        #expand("{{outdir}}/variants_{{variants}}/vireo/nclones{nclones}/{f_out}",
                    #zip,
                    #nclones=nclonelist, f_out=["enrichment/volcano_Fisher_foldNorm.png", "variants.ipynb"]),
        expand("{{outdir}}/variants_{{variants}}/vireo/nclones{nclones}/enrichment/volcano_Fisher_foldNorm.png",
                    nclones=nclonelist)
        #expand("{{outdir}}/variants_{{variants}}/vireo/nclones{nclones}/
    output:
        temp("{outdir}/clones/variants_{variants}/vireo/temp/.tmp")

    shell: "touch {output}"

################################################
## Till here
################################################


################################################
#### COMMUNITY-BASED DETECTION ALGORITHM
## So far only for mgatkdonor, not for simple
################################################
rule knn:
    input:
        "{outdir}/clones/variants_mgatkdonor/donor{d}/mgatk_donor/d{d}.variant.rds"
    output:
        "{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
        report("{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/donor{d}.variants.labels.png",
               category="enrichment")
    params:
        mgatk_in = lambda wildcards, input: input[0].replace(".variant.rds", ".af.tsv"),
        name = lambda wildcards: f"donor{wildcards.d}", #"/data2/mito_lineage/data/processed/mttrace/TcellDupi_may17_2021/MTblacklist/pre/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/pre.variant.rds"
        donor = lambda wildcards: wildcards.d,
        outdir = lambda wildcards, output: dirname(output[0]),
        kparam = lambda wildcards: wildcards.kparam,
        note = join(ROOT_DIR, "R_scripts", "knn_clones.ipynb"),
    shell:
        "papermill -k ir -p mgatk_in {params.mgatk_in} -p name {params.name} -p donor {params.donor} -p outdir {params.outdir} -p kparam {params.kparam} {params.note} {params.outdir}/output.ipynb"

rule knn_mgatkdonor_concat:
    input:
        expand("{{outdir}}/variants_mgatkdonor/knn/kparam_{{kparam}}/donor{d}/cells_meta.tsv",
               d = np.arange(config["N_DONORS"])),
    output:
         "{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/concat/cells_meta.tsv"
    run:
        all = []
        cols = ["donor", "lineage"]
        for i in input:
            print(i)
            all.append(pd.read_csv(i, sep='\t'))
            print(all[-1].head())
            if "donor_index" in all[-1].columns.values:
                cols.append("donor_index")
            if "lineage_index" in all[-1].columns.values:
                cols.append("lineage_index")
        all = pd.concat(all, ignore_index=True).sort_values(cols)
        if 'level_0' in all.columns:
                all = all.drop('level_0', axis=1)
        all.to_csv(output[0], sep='\t', index=False)


rule knn_enrichment:
    input:
        "{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/concat/cells_meta.tsv",
    output:
        report("{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/enrichment/volcano_Fisher_foldNorm.png", category="enrichment"),
    params:
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        script = join("src", "lineage", "lineage_enrichment.py"),
        samples=",".join(samples.index)
    shell: "python {params.script} {input} {params.OUTDIR} {params.samples}"


rule knn_process:
    input:
        expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{kparam}/enrichment/volcano_Fisher_foldNorm.png",
               kparam=config["params_clones"]["knn"]["params"]["resolution"])
    output: "{outdir}/clones/variants_mgatkdonor/knn/temp/.tmp"#.pipeline"
    shell:
        "touch {output}"
    #outdir <- ""#"/data2/mito_lineage/Analysis/annotation/output/data/TcellDupi_may17_2021/MTblacklist/"
        # Cluster parameters
        #
################################################

########################################################################
## Break up variants by clones and get the types of variants
rule clones_type_variants:
    input:
        "{{outdir}}/variants_{{variants}}/vireo/nclones{{nclones}}/donor{d}.labels.png"
    output: "{outdir}/clones/variants_{variants}/{method}/nclones{nclones}/variants.ipynb"
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


# Sort of like an all rule to bring downstream results together.
rule complete_lineage:
    input:
        "{outdir}/clones/variants_{variants}/{method}/temp/.tmp",
    output: "{outdir}/clones/variants_{variants}/{method}/.completed"
    shell: "touch {output}"


