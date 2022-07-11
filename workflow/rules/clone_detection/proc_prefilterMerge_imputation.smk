#configfile:  "parameters/pipeline/cosmo_server/jan21_2021.yaml"
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

params_ft = config["filters"]
ft = params_ft["params"]
samples = config["samples"]

mt_ref_fa = config["genome_path"][config["genome"]]["mt_ref_fa"]

module mgatkMod:
    snakefile: "../mgatk.smk"
    config: config

wildcard_constraints:
    sample = "^((?!merged).)*$"


rule merge_pileup_conditions:
    input:
        cov_dirs = expand("{{out_dir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/{sample}.coverage.strands.txt", sample=samples.index),
        mult_cells = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb",
    output:
        merge_pile = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/coverage_merged/donor{d}/merge.coverage.strands.txt",
                            d=np.arange(config["N_DONORS"]))
    params:
        #indir = lambda wildcards, input: dirname(input.cov_dirs[0]),
        #outdir = lambda wildcards, output: dirname(dirname(output[0])),
        samples = samples.index
    run:
        print('mult_cells', input.mult_cells)
        print(type(input.mult_cells))
        try:
            cells_meta = pd.read_csv(join(dirname(input.mult_cells),"cells_meta.tsv"), sep="\t")
        except:
            cells_meta = pd.read_csv(join(dirname(input.mult_cells[0]),"cells_meta.tsv"), sep="\t")

        from collections import defaultdict
        donor_cov_d = defaultdict(list)
        samples = params.samples
        for d, cells_don_df in cells_meta.groupby("donor"):
            for n in ["A", "G", "C", "T", "coverage"]:
                print('n', n)
                curr_pile = []
                for ind, curr_sample in enumerate(samples):
                    print(curr_sample)
                    cov_df = pd.read_csv(join(dirname(input.cov_dirs[ind]), f"{curr_sample}.{n}.strands.txt"), header=None)
                    cov_df[1] = cov_df[1].apply(lambda x: f"{x}_{curr_sample}")
                    # filter for cells in this donor
                    curr_don_cov = cov_df.loc[cov_df[1].isin(cells_don_df["ID"])].copy()
                    donor_cov_d[(d,n)].append(curr_don_cov)
                #curr_pile.append(cov_df)
        for k in donor_cov_d:
            d,n = k[0], k[1]
            print('d,n', d, n)
            curr_df = pd.concat(donor_cov_d[k], axis=0, ignore_index=True)
            print(curr_df.head())
            (curr_df).to_csv(join(dirname(output.merge_pile[int(d)]), f"merge.{n}.strands.txt"), header=False, index=False)
            donor_cov_d[k] = None


rule create_filters_prefilterMerge:
    input:
        concat_dir = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/coverage_merged/donor{d}/merge.coverage.strands.txt",
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/coverage_merged/donor{d}/out.ipynb",
        cov = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/coverage_merged/donor{d}/merged.coverage.txt",
        af = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/coverage_merged/donor{d}/af_by_cell.tsv",
        fig = report("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/coverage_merged/donor{d}/heatmap.png",
                      category="mtpreproc", subcategory="filter"),
        fig2 = report("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/coverage_merged/donor{d}/initial_cell_depth.png",
                      category="mtpreproc",subcategory="filter"),
    params:
        concat_d = lambda wildcards, input: dirname(input.concat_dir),
        ref_fa = mt_ref_fa, # config["genome_path"][config['genome']]['mt_ref_fa'],
        name = "merged", # lambda wildcards: wildcards.sample,
        #filt_params = get_filt,
        note = join(ROOT_DIR, "workflow/notebooks/af_filter/filter_af_by_cell_v02.ipynb")
    version: "v02"
    resources:
        mem_mb=90000
    shell:
        """
        papermill -p scpileup_dir {params.concat_d} -p af_f {output.af} -p mt_ref_fa {params.ref_fa} -p name {params.name} \
        -p min_cells {wildcards.mincells} -p min_reads {wildcards.minreads} -p topn {wildcards.topN} \
        -p het_thresh {wildcards.hetthresh} -p min_het_cells {wildcards.minhetcells} -p het_count_thresh {wildcards.hetcountthresh} \
        -p bq_thresh {wildcards.bqthresh} {params.note} {output.note}
        """


# rule get_refAllele_prefiltermerge:
#     #input: config["mt_ref_fa"],
#     params: params['mgatk']["chrM_refAllele"]
#     output: "{out_dir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/donor{d}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/chrM_refAllele.txt",
#     shell: 'cp {params} {output}'
#

# use rule get_refAllele from mgatkMod as preftilermerge_get_refAllele with:
#     output:
#         refAllele = "{out_dir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/donor{d}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/chrM_refAllele.txt",

rule prefiltermerge_mgatk:
    """ Run both toSeurat and call variants in one script"""
    input:
        all = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/coverage_merged/donor{d}/merged.coverage.txt",
        refAllele = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/coverage_merged/donor{d}/chrM_refAllele.txt",
    output:
        vars_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/coverage_merged/donor{d}/merged.variant.rds",
        vars_qc = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/coverage_merged/donor{d}/merged.variantQC.png",
    params:
        data_dir = lambda wildcards, input: dirname(input.all),
        sample = "merged",
        outdir = ".",
        script = join(ROOT_DIR, "R_scripts", "wrap_mgatk.R")
    shell:
        "{params.script} {params.data_dir} {params.outdir}/{params.sample} FALSE "


#rule prefiltermerge_mgatk_to_vireoIn:
use rule mgatk_to_vireoIn from mgatkMod as mergeImpute_mgatk_to_vireoIn with:
    input:
        mgatk = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/coverage_merged/donor{d}/merged.variant.rds",
    output:
        vireofmt = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/coverage_merged/donor{d}/cellSNP.tag.AD.mtx"
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        sample = "merged" # lambda wildcards: wildcards.sample
    # shell:
    #     "python src/mgatk_to_vireo.py {params.indir} {params.outdir} {params.sample}"


# use rule merge_pileup_conditions from prefilterImputeMod as imputeproc_merge_pileup_conditions with:
#     input:
#         cov_dirs = expand("{{output}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/{sample}.coverage.strands.txt", sample=samples.index),
#         mult_cells = expand("{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/multiplex/cells_meta.tsv",
#                             mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
#                             hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'])
#     output:
#         merge_pile = expand("{{output}}/data/coverage_merged/donor{d}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/merge.coverage.strands.txt",
#                             d=config["N_DONORS"])

rule scPileup_prefilterMergeImpute:
    input:
        #vars_f = "{out_dir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/merged.variant.rds",
        vireofmt = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/coverage_merged/donor{d}/cellSNP.tag.AD.mtx"
    output:
        af = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_prefilterMerge_impute/donor{d}/af.tsv",
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_prefilterMerge_impute/donor{d}/output.ipynb",
    params:
        note= join(ROOT_DIR, "workflow/notebooks/af_filter/impute_pileup_prefilterMerge.ipynb"),
        sample =  "merged", #"",".join(samples.index),
        pileup_d = lambda wildcards, input: dirname(input.vireofmt),
        outdir = lambda wildcards, output: dirname(output.af),
    shell: "papermill -p pileup_d {params.pileup_d} -p outdir {params.outdir} {params.note} {output.note}"



# module knnMod:
#     snakefile: "./rules/clone_detection/knn.smk"
#     config: config

#use rule knn from knnMod as knn_knn with:
rule knn_prefilterMergeImpute:
    input:
        af = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_prefilterMerge_impute/donor{d}/af.tsv",
    output:
        cells= "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_prefilterMerge_impute/knn/kparam_{kparam}/donors/donor{d}/cells_meta.tsv",
        fig =  "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_prefilterMerge_impute/knn/kparam_{kparam}/donors/donor{d}/donor{d}.variants.labels.png"
    params:
        indir = lambda wildcards, input: dirname(input.af),
        name = lambda wildcards: f"donor{wildcards.d}",
        outdir = lambda wildcards, output: dirname(output[0]),
        note = join(ROOT_DIR, "R_scripts", "knn_clones_init.ipynb")
    shell:
        "papermill -k ir -p indir {params.indir} -p name {params.name} -p donor {wildcards.d} -p outdir {params.outdir} -p kparam {wildcards.kparam} {params.note} {params.outdir}/output.ipynb"


#use rule knn_concat from knnMod as knn_knn_concat with:
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


rule knn_concat_prefilterMergeImpute:
    input:
        cells = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_prefilterMerge_impute/knn/kparam_{{kparam}}/donors/donor{d}/cells_meta.tsv",
                       d=np.arange(config["N_DONORS"])),
    output:
        cells_all = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_prefilterMerge_impute/knn/kparam_{kparam}/cells_meta.tsv"
    run: concat_knn(input, output[0])

