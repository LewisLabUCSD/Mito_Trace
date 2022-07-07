from src.config import ROOT_DIR
import pandas as pd
from os.path import join, dirname
import os
import numpy as np
from icecream import ic

# def get_counts_in(wildcards):
#     w = wildcards
#     if w.variants == "mgatkdonor":
#         return f"{w.output}/data/merged/MT/cellr_{w.cellrbc}/numread_{w.num_read}/filters/minC{w.mincells}_minR{w.minreads}_topN{w.topN}_hetT{w.hetthresh}_hetC{w.minhetcells}_hetCount{w.hetcountthresh}_bq{w.bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{w.d}/mgatk/d{w.d}.variant.rds"
#     elif w.variants == "simple":
#         return f"{w.output}/data/merged/MT/cellr_{w.cellrbc}/numread_{w.num_read}/filters/minC{w.mincells}_minR{w.minreads}_topN{w.topN}_hetT{w.hetthresh}_hetC{w.minhetcells}_hetCount{w.hetcountthresh}_bq{w.bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb"
#     elif w.variants == "init":
#         return f"{w.output}/data/merged/MT/cellr_{w.cellrbc}/numread_{w.num_read}/filters/minC{w.mincells}_minR{w.minreads}_topN{w.topN}_hetT{w.hetthresh}_hetC{w.minhetcells}_hetCount{w.hetcountthresh}_bq{w.bqthresh}/mgatk/vireoIn/multiplex/clones_init/donor{w.d}/af.tsv"
#     return




#samples = config["samples"] #pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)

########################################################################
## Variant Intermediate step for other workflows.
## Extract pileups for each donor from original filter matrices and run mgatk
########################################################################
# def get_coverage(wildcards):
#     return f"{config['cov_indir']['sample']}/{wildcards.sample}.coverage.txt"
rule convert_to_af:
    input:
        cells_meta = "{outdir}/cells_meta.tsv",
        counts_in = "af.tsv"
    output:
        note="{outdir}/sc_af/donor{d}/sc_af.ipynb"
        #af="{outdir}/sc_af/donor{d}/af.tsv",
    params:
        note = join(ROOT_DIR, "workflow", "notebooks", "clone_af_dendrograms", "Convert_to_AF.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
        indir = lambda wildcards, input: dirname(input.cells_meta),
        counts_indir = lambda wildcards, input: dirname(input.counts_in),
        var_type = "init"
    shell:
        "papermill -p INDIR {params.indir} -p COUNT_INDIR {params.counts_indir} -p OUTDIR {params.outdir} -p DONOR {wildcards.d}  -p var_type {params.var_type} {params.note} {output.note}"


rule barcodes_btwnClones:
    input:
        cells_meta = "{outdir}/cells_meta.tsv",
        af_note = "{outdir}/sc_af/donor{d}/sc_af.ipynb",
    output:
        note = "{outdir}/barcodes/btwnClones/donor{d}.ipynb",
        # dendros = report("{outdir}/barcodes/btwnClones/donor{d}.na.clust.max2.AF.png",
        #                   category="lineage", subcategory="Clone Barcodes (conditions separate): Donor {d}"),
        # dendrodp  =  report("{outdir}/barcodes/btwnClones/donor{d}.na.clust.max2.DP.png",
        #     category="lineage", subcategory="Clone Barcodes (conditions separate): Donor {d}"),
        # dendrosCond = report("{outdir}/barcodes/btwnClones/donor{d}.NoCondition.na.clust.max2.AF.png",
        #     category="lineage",subcategory="Clone Barcodes: Donor {d}"),
        # dendrodpCond =  report("{outdir}/barcodes/btwnClones/donor{d}.NoCondition.na.clust.max2.DP.png",
        #     category="lineage", subcategory="Clone Barcodes: Donor {d}"),
    params:
        note = join(ROOT_DIR, "workflow", "notebooks", "clone_af_dendrograms", "MT_btwnClones_Barcode.ipynb"),
        indir = lambda wildcards, input: dirname(input.cells_meta),
        outdir = lambda wildcards, output: dirname(output.note),
    shell: "papermill -p INDIR {params.indir} -p OUTDIR {params.outdir} -p DONOR {wildcards.d}  {params.note} {output.note}"



rule barcodes_btwnClones_dendro:
    input:
        cells_meta = "{outdir}/cells_meta.tsv",
        af_note = "{outdir}/sc_af/donor{d}/sc_af.ipynb",
    output:
        note = "{outdir}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.ipynb",
        mean = "{outdir}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.mean.csv",
        res = report(multiext("{outdir}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.",
                                  "clones_dendro.csv", "dendrogram_pvals.txt",
                                  "dendro.NoCondition.max2.AF.png"),
                         category="lineage")

    params:
        note = join(ROOT_DIR, "workflow", "notebooks", "clone_af_dendrograms", "MT_btwnClones_Barcode_dendro.ipynb"),
        indir = lambda wildcards, input: dirname(input.cells_meta),
        outdir = lambda wildcards, output: dirname(output.note),
        dendro_thresh = 0.6
    shell: "papermill -p INDIR {params.indir} -p OUTDIR {params.outdir} -p DONOR {wildcards.d} -p dendroThresh {wildcards.dendro_thresh} {params.note} {output.note}"



rule barcodes_inClones:
    input:
        cells_meta = "{outdir}/cells_meta.tsv",
        af_note = "{outdir}/sc_af/donor{d}/sc_af.ipynb",
    output:
        note = "{outdir}/barcodes/inClones/donor{d}.ipynb",
    params:
        note = join(ROOT_DIR, "workflow/notebooks/clone_af_dendrograms", "MT_inClones_Barcode.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
        indir = lambda wildcards, input: dirname(input.cells_meta),
    shell: "papermill -p INDIR {params.indir} -p OUTDIR {params.outdir} -p DONOR {wildcards.d} {params.note} {output.note}"


rule distinguishing_vars:
    input:
        cells_meta = "{outdir}/cells_meta.tsv",
        af_note = "{outdir}/sc_af/donor{d}/sc_af.ipynb",
    output:
        note = "{outdir}/distinct_variants/donor{d}/output.ipynb",
    params:
        note = join(ROOT_DIR, "workflow/notebooks/clone_vars", "distinguishing_vars.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
        indir = lambda wildcards, input: dirname(input.cells_meta),
    shell: "papermill -p INDIR {params.indir} -p OUTDIR {params.outdir} -p DONOR {wildcards.d} {params.note} {output.note}"


rule finalize:
    input:
        #inClones =  expand("{{outdir}}/barcodes/btwnClones/donor{d}.ipynb", d=np.arange(config["N_DONORS"])),
        btwnClones = expand("{{outdir}}/barcodes/inClones/donor{d}.ipynb", d=np.arange(config["N_DONORS"])),
        distinct_vars = expand("{{outdir}}/distinct_variants/donor{d}/output.ipynb", d=np.arange(config["N_DONORS"])),
        dendroClones = expand("{{outdir}}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.mean.csv",
                              d=np.arange(config["N_DONORS"]), dendro_thresh=[config["dendro_thresh"]]),
    output:
        out = "{outdir}/barcodes/_clone_complete.txt",
    shell: "touch {output.out}"


################################################
#### COMMUNITY-BASED DETECTION ALGORITHM
## So far only for mgatkdonor, not for simple
################################################
# rule knn:
#     input:
#         in_vars="{outdir}/clones/variants_mgatkdonor/donor{d}/mgatk_donor/d{d}.variant.rds"
#     output:
#         cells="{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
#         fig=report("{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/donor{d}.variants.labels.png",
#                category="enrichment")
#     params:
#         mgatk_in = lambda wildcards, input: input[0].replace(".variant.rds", ".af.tsv"),
#         name = lambda wildcards: f"donor{wildcards.d}", #"/data2/mito_lineage/data/processed/mttrace/TcellDupi_may17_2021/MTblacklist/pre/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/pre.variant.rds"
#         donor = lambda wildcards: wildcards.d,
#         outdir = lambda wildcards, output: dirname(output[0]),
#         kparam = lambda wildcards: wildcards.kparam,
#         note = join(ROOT_DIR, "R_scripts", "knn_clones.ipynb"),
#     shell:
#         "papermill -k ir -p mgatk_in {params.mgatk_in} -p name {params.name} -p donor {params.donor} -p outdir {params.outdir} -p kparam {params.kparam} {params.note} {params.outdir}/output.ipynb"
#
# rule knn_mgatkdonor_concat:
#     input:
#         cells=expand("{{outdir}}/clones/variants_mgatkdonor/knn/kparam_{{kparam}}/donor{d}/cells_meta.tsv",
#                d = np.arange(config["N_DONORS"])),
#     output:
#          cells_all="{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/concat/cells_meta.tsv"
#     run:
#         all = []
#         cols = ["donor", "lineage"]
#         for i in input:
#             print(i)
#             all.append(pd.read_csv(i, sep='\t'))
#             print(all[-1].head())
#             if "donor_index" in all[-1].columns.values:
#                 cols.append("donor_index")
#             if "lineage_index" in all[-1].columns.values:
#                 cols.append("lineage_index")
#         all = pd.concat(all, ignore_index=True).sort_values(cols)
#         if 'level_0' in all.columns:
#                 all = all.drop('level_0', axis=1)
#         all.to_csv(output[0], sep='\t', index=False)
#
#
# rule knn_enrichment:
#     input:
#         "{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/concat/cells_meta.tsv",
#     output:
#         report("{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/enrichment/volcano_Fisher_foldNorm.png", category="enrichment"),
#     params:
#         OUTDIR = lambda wildcards, output: dirname(output[0]),
#         script = join("src", "lineage", "lineage_enrichment.py"),
#         samples=",".join(samples.index)
#     shell: "python {params.script} {input} {params.OUTDIR} {params.samples}"
#
#
# rule knn_process:
#     input:
#         expand("{{outdir}}/clones/variants_mgatkdonor/knn/kparam_{kparam}/enrichment/volcano_Fisher_foldNorm.png",
#                kparam=config["params_clones"]["knn"]["params"]["resolution"])
#     output: "{outdir}/clones/variants_mgatkdonor/knn/temp/.tmp"#.pipeline"
#     shell:
#         "touch {output}"
#     #outdir <- ""#"/data2/mito_lineage/Analysis/annotation/output/data/TcellDupi_may17_2021/MTblacklist/"
#         # Cluster parameters
#         #
#
# rule knn_process_simple:
#     """ dummy rule to not use simple for knn
#     """
#     # input:
#     #     expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_simple/knn/kparam_{kparam}/enrichment/volcano_Fisher_foldNorm.png",
#     #            kparam=params_clones["knn"]["params"]["resolution"])
#     output: "{outdir}/clones/variants_simple/knn/temp/.tmp"#.pipeline"
#     shell:
#         "touch {output}"

################################################

########################################################################
## Break up variants by clones and get the types of variants
# rule clones_type_variants:
#     input:
#         "{{outdir}}/clones/variants_{{variants}}/vireo/nclones{{nclones}}/donor{d}.labels.png"
#     output: "{outdir}/clones/variants_{variants}/{method}/nclones{nclones}/variants.ipynb"
#     params:
#         notebook=join("src", "vireo", join("6_MT_Clones_variantTypes.ipynb")),
#         INDIR = lambda wildcards, input: dirname(input[0]),
#         OUTDIR = lambda wildcards, output: dirname(output[0]),
#         N_DONORS = config["N_DONORS"],
#         sample_names = ','.join(samples.index), # make it as a list
#         nclones = lambda  wildcards: wildcards.nclones, #lambda wildcards: wildcards.nclones, #config['multiplex']["n_clone_list"],#lambda wildcards: wildcards.nclones,
#         var_thresh=0.001,
#         vars_to_plot=10
#     shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR}  -p nclones {params.nclones} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} -p var_thresh {params.var_thresh} -p vars_to_plot {params.vars_to_plot} {params.notebook} {output} && jupyter nbconvert --to pdf {output}"
#

# Sort of like an all rule to bring downstream results together.
# rule complete_lineage:
#     input:
#         "{outdir}/clones/variants_{variants}/{method}/temp/.tmp",
#     output: "{outdir}/clones/variants_{variants}/{method}/.completed"
#     shell: "touch {output}"
#
#
