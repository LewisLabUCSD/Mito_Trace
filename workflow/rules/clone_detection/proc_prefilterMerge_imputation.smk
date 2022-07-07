""" Workflow for merge before filter and data imputation

Brief overview: Create donor folders and merge across all conditions,
using the pre-filter coverage. Then we run the filter and MGATK on each donor to
get the variants to use. This should have the imputed data, and from
that we get an af.tsv and dp.tsv files, and we run clone detection on that.
After clone detection we output the sc_af.tsv for the cells with clone
assignment (similar to before but some may be filtered). With this,
we can continue with the rest of the pipeline using the 'clones' folder
and the 'multiplex/variants_{impute_preMerge}' folder.


Workflow:
V01:
1. merge_donors_prefilter
Input:
    a. cells_meta from donor assignment
    b. coverage files from temp("{outdir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/{sample}.coverage.strands.txt") in mt_preprocess.smk
Output:
   Folder for each donor d with:
   a. cells coverage similar to b in input (i.e. input for create_filters_v02 function). Each cell has suffix with condition
2. Filter
3. MGATK and mgatk_to_vireoIn
4. Convert MGATK and coverage to af.tsv for each donor:
    Need to create it, output similar to scPileup_simpleUnion_v02
5. Clone detect:
Input:
    a) multiplex/clones_init/donor{d}/af.tsv
Output:
    cells = clones/variants_simpleUnion/knn/kparam_{{kparam}}/donor{d}/cells_meta.tsv",
6. Get sc_af for clones
7. Run everything else as is with the multiplex and clones folders.


"""

#from src.config import ROOT_DIR
#workdir: ROOT_DIR
#wildcard_constraints:
    # cellr='True|False',
    # kparam='[0-9]+',  #'%d', #[0-9]+',
    # variants = "simple|mgatkdonor|init|simpleUnion|prefilterMerge_impute",
    # d = "[0-9]+"


from src.config import ROOT_DIR
from src.utils.parse_config import read_config_file
import os
import numpy as np
from os.path import join, dirname
import pandas as pd
from snakemake.utils import min_version
min_version("6.0")
print('config', config)
import pickle

#res = join(config["outdir"], "pipeline", config["prefix"])

samples =  config["samples"] #pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)
params = config
#anno_res = join(config["outdir"], "annotation", "data", params['annotations']['version'], config["prefix"])

################################################################
## Import from prior snakefile modules
## Here, we redefine the input to be based on our config['files'] dictionary
from snakemake.utils import min_version
min_version("6.0")


# module knnMod:
#     snakefile: "./knn.smk"
#     config: config
#

#use rule * from mtpreprocMod as mtpreproctwo_*
#use rule * from mgatkMod as mgatktwo_*
#use rule * from multMod as multtwo_*

#use rule * from knnMod as knntwo_*

print('samples', samples.index)
params_ft = config["filters"]
ft = params_ft["params"]
# rule merge_pileup_conditions:
#     input:
#         cov_dirs = expand("{{outdir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/{sample}.coverage.strands.txt", sample=samples.index),
#         mult_cells = expand("{{outdir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/multiplex/cells_meta.tsv",
#                             mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
#                             hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'])
#     output: merge_pile = expand("{{outdir}}/data/coverage_merged/donor{d}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/merge.coverage.strands.txt", d=config["N_DONORS"])
#     params:
#         indir = lambda wildcards, input: dirname(input.cov_dirs[0]),
#         outdir = lambda wildcards, output: dirname(dirname(output[0])),
#         samples = samples.index
#     run:
#         cells_meta = pd.read_csv(input.mult_cells, sep="\t")
#         from collections import defaultdict
#         donor_cov_d = defaultdict(list)
#         for d, cells_don_df in cells_meta.groupby("donor"):
#             for n in ["A", "G", "C", "T", "coverage"]:
#                 print('n', n)
#                 curr_pile = []
#                 for ind, curr_sample in enumerate(samples):
#                     cov_df = pd.read_csv(join(params.indir, f"{curr_sample}.{n}.strands.txt"), header=None)
#                     cov_df[1] = cov_df[1].apply(lambda x: f"{x}_{curr_sample}")
#                     # filter for cells in this donor
#                     curr_don_cov = cov_df.loc[cov_df[1].isin(cells_don_df["raw ID"])].copy()
#                     donor_cov_d[(d,n)].append(curr_don_cov)
#                 #curr_pile.append(cov_df)
#         for (d,n), val in donor_cov_d:
#             (pd.concat(val, axis=0)).to_csv(join(params.outdir, f"donor{d}", f"merge.{n}.strands.txt"), sep="\t")

module mtpreprocMod:
    snakefile: "../mt_preprocess.smk"
    config: config

module mgatkMod:
    snakefile: "../mgatk.smk"
    config: config

module multMod:
    snakefile: "../multiplex.smk"
    config: config


# use rule create_filters_v02 from mtpreprocMod as mtpreproctwo_create_filters_v02 with:
#     input:
#         concat_dir = "{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/merge.coverage.strands.txt"
#     output:
#         note = "{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/out.ipynb",
#         cov = "{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/merged.coverage.txt",
#         af = "{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/af_by_cell.tsv",
#         fig = report("{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/heatmap.png",
#                       category="mtpreproc", subcategory="filter"),
#         fig2 = report("{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/initial_cell_depth.png",
#                       category="mtpreproc",subcategory="filter"),
#     params:
#         ref_fa = config["genome_path"][config['genome']]['mt_ref_fa'],
#         name = "merge",
#         concat_d = lambda wildcards, input: dirname(input.concat_dir),
#         note = join(ROOT_DIR, "workflow/notebooks/af_filter/filter_af_by_cell_v02.ipynb")
#

# rule create_filters_prefilterMerge:
#     input:
#         concat_dir = "{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/merge.coverage.strands.txt"
#     output:
#         note = "{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/out.ipynb",
#         cov = "{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/merged.coverage.txt",
#         af = "{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/af_by_cell.tsv",
#         fig = report("{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/heatmap.png",
#                       category="mtpreproc", subcategory="filter"),
#         fig2 = report("{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/initial_cell_depth.png",
#                       category="mtpreproc",subcategory="filter"),
#     params:
#         concat_d = lambda wildcards, input: dirname(input.concat_dir),
#         ref_fa = config["genome_path"][config['genome']]['mt_ref_fa'],
#         name = lambda wildcards: wildcards.sample,
#         #filt_params = get_filt,
#         note = join(ROOT_DIR, "workflow/notebooks/af_filter/filter_af_by_cell_v02.ipynb")
#     version: "v02"
#     resources:
#         mem_mb=90000
#     shell:
#         """
#         papermill -p scpileup_dir {params.concat_d} -p af_f {output.af} -p mt_ref_fa {params.ref_fa} -p name {params.name} \
#         -p min_cells {wildcards.mincells} -p min_reads {wildcards.minreads} -p topn {wildcards.topN} \
#         -p het_thresh {wildcards.hetthresh} -p min_het_cells {wildcards.minhetcells} -p het_count_thresh {wildcards.hetcountthresh} \
#         -p bq_thresh {wildcards.bqthresh} {params.note} {output.note}
#         """

#
# use rule mgatk from mgatkMod as mgatktwo_mgatk with:
#     input:
#         cov = "{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/donor{d}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/merged.coverage.txt",
#         refAllele = "{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/donor{d}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/chrM_refAllele.txt",
#     output:
#         vars_f = "{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/donor{d}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/merged.variant.rds",
#         vars_qc = "{outdir}/data/coverage_merged/donor{d}/  MT/cellr_{cellrbc}/numread_{num_read}/donor{d}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/merged.variantQC.png",


########################################################################
## bring back to initial folders in the output here, changing the
## folder from {outdir}/data/coverage_merged to {outdir}/data/merged.
########################################################################
# create multiplex/variants_{prefilterMerge_impute}/donor{d}/af.tsv & dp.tsv Convert MGATKT output
# impute here
# rule mgatk_to_af:
#     input:
#         vars_f = "{outdir}/data/coverage_merged/donor{d}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/merged.variant.rds",
#     output:
#         af = "{outdir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_prefilterMerge_impute/donor{d}/af.tsv",
#         note = "{outdir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_prefilterMerge_impute/donor{d}/output.ipynb",
#     params:
#         note= join(ROOT_DIR, "workflow/notebooks/af_filter/impute_pileup.ipynb")
#     shell: "papermill {params.note} {output.note}"


# use rule knn from knnMod as knntwo_knn with:
#     input:
#         af = "{outdir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_prefilterMerge_impute/donor{d}/af.tsv",
#     output:
#         cells= "{outdir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_prefilterMerge_impute/knn/kparam_{kparam}/donors/donor{d}/cells_meta.tsv",
#         fig =  "{outdir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_prefilterMerge_impute/knn/kparam_{kparam}/donors/donor{d}/donor{d}.variants.labels.png"
#
#
#
# use rule knn_concat from knnMod as knntwo_knn_concat with:
#     input:
#         cells = expand("{{outdir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_prefilterMerge_impute/knn/kparam_{{kparam}}/donors/donor{d}/cells_meta.tsv",
#                        d=np.arange(config["N_DONORS"])),
#     output:
#         cells_all = "{outdir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_prefilterMerge_impute/knn/kparam_{kparam}/cells_meta.tsv"
#
