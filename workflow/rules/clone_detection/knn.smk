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

########################################################################
## Variant Intermediate step for other workflows.
## Extract pileups for each donor from original filter matrices and run mgatk
########################################################################
# def get_coverage(wildcards):
#     return f"{config['cov_indir']['sample']}/{wildcards.sample}.coverage.txt"
def concat_knn(in_files, out_name):
    all = []
    cols = ["donor", "lineage"]
    for i in in_files:
        #print(i)
        all.append(pd.read_csv(i,sep='\t'))
        #print(all[-1].head())
        if "donor_index" in all[-1].columns.values:
            cols.append("donor_index")
        if "lineage_index" in all[-1].columns.values:
            cols.append("lineage_index")
    all = pd.concat(all,ignore_index=True).sort_values(cols)
    if 'level_0' in all.columns:
        all = all.drop('level_0',axis=1)
    all.to_csv(out_name,sep='\t',index=False)


def get_knn_script(cfg):
    if cfg["variants_type"] == "init":
        join(ROOT_DIR, "R_scripts", "knn_clones_init.ipynb"),
    elif cfg["variants_type"] == "mgatkdonor":
        return join(ROOT_DIR, "R_scripts", "knn_clones.ipynb"),

rule knn:
    input:
        af ="{outdir}/donor{d}/af.tsv"
    output:
        cells= "{outdir}/knn/kparam_{kparam}/donors/donor{d}/cells_meta.tsv",
        fig = "{outdir}/knn/kparam_{kparam}/donors/donor{d}/donor{d}.variants.labels.png"
    params:
        indir = lambda wildcards, input: dirname(input.af),
        name = lambda wildcards: f"donor{wildcards.d}",
        outdir = lambda wildcards, output: dirname(output[0]),
        note = get_knn_script(config) #join(ROOT_DIR, "R_scripts", "knn_clones_init.ipynb"),
    shell:
        "papermill -k ir -p indir {params.indir} -p name {params.name} -p donor {wildcards.d} -p outdir {params.outdir} -p kparam {wildcards.kparam} {params.note} {params.outdir}/output.ipynb"


rule knn_concat:
    input:
        cells= expand("{{outdir}}/knn/kparam_{{kparam}}/donors/donor{d}/cells_meta.tsv",
                       d=np.arange(config["N_DONORS"])),
    output:
        cells_all="{outdir}/knn/kparam_{kparam}/cells_meta.tsv"
    run: concat_knn(input, output[0])


#
# rule knn_mgatkdonor:
#     input:
#         in_vars="{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/mgatk/d{d}.variant.rds"
#     output:
#         cells= "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
#         fig = "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/donor{d}.variants.labels.png"
#     # fig=report("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/donor{d}.variants.labels.png",
#     #        category="enrichment")
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
#
#
# rule knn_mgatkdonor_concat:
#     input:
#         cells=expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{{kparam}}/donor{d}/cells_meta.tsv",
#             d = np.arange(config["N_DONORS"])),
#     output:
#         cells_all="{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{kparam}/concat/cells_meta.tsv"
#     run:
#         all = []
#         cols = ["donor", "lineage"]
#         for i in input:
#             #print(i)
#             all.append(pd.read_csv(i, sep='\t'))
#             #print(all[-1].head())
#             if "donor_index" in all[-1].columns.values:
#                 cols.append("donor_index")
#             if "lineage_index" in all[-1].columns.values:
#                 cols.append("lineage_index")
#         all = pd.concat(all, ignore_index=True).sort_values(cols)
#         if 'level_0' in all.columns:
#             all = all.drop('level_0', axis=1)
#         all.to_csv(output[0], sep='\t', index=False)
#
# rule knn_init:
#     input:
#         af ="{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_init/donor{d}/af.tsv"
#     output:
#         cells= "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_init/knn/kparam_{kparam}/donors/donor{d}/cells_meta.tsv",
#         fig = "{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_init/knn/kparam_{kparam}/donors/donor{d}/donor{d}.variants.labels.png"
#     # fig=report("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{kparam}/donors/donor{d}/donor{d}.variants.labels.png",
#     #        category="enrichment")
#     params:
#         indir = lambda wildcards, input: dirname(input[0]),
#         name = lambda wildcards: f"donor{wildcards.d}", #"/data2/mito_lineage/data/processed/mttrace/TcellDupi_may17_2021/MTblacklist/pre/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/pre.variant.rds"
#         donor = lambda wildcards: wildcards.d,
#         outdir = lambda wildcards, output: dirname(output[0]),
#         kparam = lambda wildcards: wildcards.kparam,
#         note = join(ROOT_DIR, "R_scripts", "knn_clones_init.ipynb"),
#     shell:
#         "papermill -k ir -p indir {params.indir} -p name {params.name} -p donor {params.donor} -p outdir {params.outdir} -p kparam {params.kparam} {params.note} {params.outdir}/output.ipynb"
#
#
# rule knn_init_concat:
#     input:
#         cells=expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_init/knn/kparam_{{kparam}}/donors/donor{d}/cells_meta.tsv",
#             d=np.arange(config["N_DONORS"])),
#     output:
#         cells_all="{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_init/knn/kparam_{kparam}/cells_meta.tsv"
#     run: concat_knn(input, output[0])


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
rule complete_lineage:
    input:
        "{outdir}/clones/variants_{variants}/{method}/temp/.tmp",
    output: "{outdir}/clones/variants_{variants}/{method}/.completed"
    shell: "touch {output}"


