rule rclone:
    input:
         "{input}"
    output:
        "{input}/output/rsync_copy.ok"
    params:
        rsync_loc = config["rsync_location"]
    shell: "rcopy {input} drive:{rsync_loc} > rsync_copy.ok && mv rsync_copy.ok {output}/"

report: "report/workflow.rst"

#from src.config import ROOT_DIR
#workdir: ROOT_DIR
wildcard_constraints:
    cellr='True|False',
    kparam='[0-9]+',  #'%d', #[0-9]+',
    variants = "simple|mgatkdonor|init|simpleUnion",
    d = "[0-9]+"

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

########################################################################
# Setup parameters and outdir
########################################################################
res = join(config["outdir"], "pipeline", config["prefix"])

params = read_config_file(config["config"])
samples = pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)
#anno_res = join(config["outdir"], "annotation", "data", params['annotations']['version'], config["prefix"])
anno_res = join(config["outdir"], "annotation", "data", config["prefix"])


def create_clone_report(inputs, output):
    return

###################################
## Rules and workflow
###################################
###################################
rule all:
    input:
        expand("{out_dir}/reports/big_data_{version}/clones/clones.pdf",
            out_dir=res, version=params["version"]),
        expand("{out_dir}/reports/big_data_{version}/_complete.txt", out_dir=res, version=params["version"])


rule clone_report:
    input:
        anno_dir = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/",
        clone_dir = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/",
        annClo_dir = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/",
        aggreg_dir = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/merge_enrich_and_lineage",
        hypergeom_dir = "{out_dir}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/",
    output:
        "{out_dir}/reports/big_data_{version}/clones/clones.pdf",
        note = "{out_dir}/reports/big_data_{version}/clones/clones.pdf",
    params:
        rscript = "workflow/notebooks/summarize_clones.ipynb",
        params_f = config["params"]
    run:
        "papermill -p {params.rscript} {output.note}"
        #create_clone_report(input, output)

#     input:
#         # A. Preproc: Initial MT coverage per position per cell
#         expand("{out_dir}/figures/{sample}/MT/cellr_{cellrbc}/numread_{num_read}_MT_position_coverage.png",
#             out_dir=res,sample=samples.index, num_read=num_reads_filter, cellrbc=cellrbc),
#
#         # B. Multiplex: Multiplexed VAF and depth dendrograms for each donor and selected variants
#         expand("{out_dir}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/multiplex.pdf",
#             out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
#             mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
#             hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh']),
#         expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_init/donor{d}/af.tsv",
#             out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
#             mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
#             hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
#             method=params_clones["method"],
#             kparam=params_clones["knn"]["params"]["resolution"],d=np.arange(config["N_DONORS"])),
#         # Variants
#         expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_variants/anno_variants.tsv",
#             out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
#             mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
#             hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh']),
#
#         # Clones (lineage)
#         # donor specific
#         # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
#         #     out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
#         #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
#         #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'], variants=[x for x in params_clones["variants"] if x!="simple"],
#         #     kparam=params_clones["knn"]["params"]["resolution"],d=np.arange(config["N_DONORS"])),
#         # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_init/knn/kparam_{kparam}/cells_meta.tsv",
#         #     out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
#         #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
#         #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
#         #     method=params_clones["method"],
#         #     kparam=params_clones["knn"]["params"]["resolution"],
#             #variants=[x for x in params_clones["variants"] if x!="simple"]
#        # ),
#         # Clone stats
#         expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/_clone_complete.txt",
#             out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
#             mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
#             hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
#             method=params_clones["method"],
#             kparam=params_clones["knn"]["params"]["resolution"],
#             variants=[x for x in params_clones["variants"] if x!="simple"]
#         ),
#
#         # # Enrichment
#         expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/_enrichment_complete",
#             out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
#             variants=[x for x in params_clones["variants"] if x != "simple"],
#             mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
#             hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
#             kparam=params_clones["knn"]["params"]["resolution"]),
#

#         # MT Plots
#         expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/thr__{t}_rt__{rt}/annotation_clones/de_clone_btwnvars_RNA_af/mtPlots.ipynb",
#             out_dir=res,cellrbc=cellrbc,num_read=num_reads_filter,
#             mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
#             hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
#             kparam=config["params_clones"]["knn"]["params"]["resolution"],
#             variants=[x for x in params_clones["variants"] if x != "simple"], gff=gff,
#             t=params_annclo["t"], rt=params_annclo["rt"]),
#
#         # Cluster DE
#
#         # Condition DE
#         ##  Methods compare
#         expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/comparisons/comparisons.ipynb",
#             out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
#             mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
#             hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh']),
#
#         # # C. Clones: Variant types
#         # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb",
#         #        output=res, cellrbc=cellrbc, num_read=num_reads_filter,
#         #        mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
#         #        hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], n_clones=config['multiplex']['n_clone_list']),

################################################################
## Import from prior snakefile modules
## Here, we redefine the input to be based on our config['files'] dictionary
from snakemake.utils import min_version
min_version("6.0")
module mtpreprocMod:
    snakefile: "./rules/mt_preprocess.smk"
    config: params

module mgatkMod:
    snakefile: "./rules/mgatk.smk"
    config: params

module multMod:
    snakefile: "./rules/multiplex.smk"
    config: params

use rule * from mtpreprocMod
use rule * from mgatkMod

########################################################################
## 4. De-multiplex the merged conditions using Vireo
########################################################################
use rule * from multMod as mult_*
# ########################################################################
# ## Variants A
# ########################################################################
module procDonMgatkMod:
    snakefile: "./rules/clone_detection/proc_donormgatk.smk"
    config: config
use rule * from procDonMgatkMod as donMGATK_*

# ########################################################################
# Variants B
# ########################################################################
module procInitMod:
    snakefile: "./rules/clone_detection/proc_init.smk"
    config: config

use rule * from procInitMod as procInitMod_*

########################################################################
## Variants C: 'Simple': Directly from multiplex output.
## Variants not seen in one condition might be 0'd out, which can bias
## the clone calling. Only in Vireo
########################################################################

########################################################################
## Variants D: 'simpleUnion': Union of variants from both mgatk outputs. Similar to init,
# but it also adds in coverage if a variant was filtered out in a specific condition.
########################################################################

module procSimUnMod:
    snakefile: "./rules/clone_detection/proc_simpleUnion.smk"
    config: config


use rule * from procSimUnMod as procSimUnMod_*

module cloneMod:
    snakefile: "./rules/clone_detection/clone.smk"
    config: config

use rule * from cloneMod as clone_*


########################################################################
## Enrichment
########################################################################
module cloneEnrichMod:
    snakefile: "./rules/clone_detection/enrichment.smk"
    config: params

use rule * from cloneEnrichMod as cloneEnrich_*


############################################
## Annotation
############################################
module annoMod:
    snakefile: "./rules/annotation_pipe.smk"
    config: config

use rule * from annoMod as anno_*

use rule integrateSignac from annoMod as anno_integrateSignac with:
    input:
       a = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/allSamples.rds",
    output:
        a = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.ipynb",
        b = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds",
        c = report(expand("{{out_dir}}/data/annotation/gff_{{gff}}/mergedSamples/{f}",
                      f=["QC_02.png",
                         "integrated.merged.compare.png",
                         "integrated.batch.png",
                         "integrated.lsi.clusters.png"]), category="Nuclear ATAC", subcategory="UMAP")


############################################
## Clones DE
############################################
module annCloMod:
    snakefile: "./rules/annotation_clones.smk"
    config: params

use rule * from annCloMod as annClo_*

use rule addClones from annCloMod as annClo_addClones with:
    input:
        noc = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds",
        clones = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/cells_meta.tsv"
    output:
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/SE.rds",
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/addClones.ipynb"

use rule addClones from annCloMod as knn_annClo_addClones with:
    input:
        noc = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds",
        clones = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv"
    output:
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/SE.rds",
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/addClones.ipynb"


## DE without clone
module btwnClustersMod:
    snakefile: "./rules/de/btwnClusters.smk"
    config: params
use rule * from btwnClustersMod as btwnClusters_*

module btwnConditionsMod:
    snakefile: "./rules/de/btwnConditions.smk"
    config: params
use rule * from btwnConditionsMod as btwnConditionsMod_*

## DE with clone
module inCloneBtwnConditionsMod:
    snakefile: "./rules/de/inClones_btwnConditions.smk"
    config: params

use rule * from inCloneBtwnConditionsMod as inCloneBtwnConditionsMod_*



use rule inVar_btwnCond_DE from annCloMod as annClo_inVar_btwnCond_DE with:
    input:
        mt_cells = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/mt_clones/thr__{t}_rt__{rt}/cells_meta.tsv",
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/SE.rds",
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/thr__{t}_rt__{rt}/annotation_clones/de_clone_btwnvars_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/de.ipynb",




############################################
############################################
rule clones_donors_table:
    input:
        cells_meta_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/nclones{nclones}/cells_meta.tsv"
    output:
        summary_table = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/nclones{nclones}/summary/counts_summary.csv",
        summary_fig = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/nclones{nclones}/summary/counts_summary.pdf"
    params:
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        notebook = join(ROOT_DIR, "src", "clones")
    shell: "papermill -p cells_meta_f {input.cells_meta_f} -p OUTDIR {params.OUTDIR} {params.notebook} {params.OUTDIR}/output.ipynb"


use rule output_largeclones from annCloMod as annClo_output_largeclones with:
    input:
        de = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb",
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/se.clonesfilt.rds",
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/TFactivity.ipynb",
        out = report(multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/",
            "dotplot.allDonors.clones.topPvals.pdf","dotplot.allDonors.clones.topPvals.noZ.pdf",
            "clone.TFactivity.topPvals.txt"),
            category="lineage and nuclear", subcategory="Method: vireo Variants: {variants} nClones: {nclones} minPct: {minPct} logThresh: {logThresh}")

## Use annotation_clones to run DE on clones
use rule output_largeclones from annCloMod as knn_annClo_output_largeclones with:
    input:
        de = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb",
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/se.clonesfilt.rds"
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/TFactivity.ipynb",
        out = report(multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/",
            "dotplot.allDonors.clones.topPvals.pdf","dotplot.allDonors.clones.topPvals.noZ.pdf",
            "clone.TFactivity.topPvals.txt"),
            category="lineage and nuclear", subcategory="Method: knn Variants: {variants} resolution: {kparam} minPct: {minPct} logThresh: {logThresh}"),


############################################
## clone expansion vs regressed DE
############################################
module cloneChangeMod:
    snakefile: "./rules/de/de_clones_change.smk"
    config: params #params_clch

use rule * from cloneChangeMod as cloneChange_*


use rule preprocess from cloneChangeMod as cloneChange_preprocess with:
    input:
        enrichment_dir = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/enrichmentNorm.csv"
    output:
        clone_change_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/clone_change.csv",
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/preproc.ipynb"


use rule btwnClone_DE from cloneChangeMod as cloneChange_btwnClone_DE with:
    input:
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/SE.rds",
        clone_change_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/clone_change.csv",
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/{de_group}/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/de.ipynb"


## get clone meta from stats and immune clusters
rule merge_enrich_and_lineage:
    input:
        clone_change_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/clone_change.csv",
        dominant_cluster = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/dominant_clone_clust/combinedConditions_dominant_cluster.csv"
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/merge_enrich_and_lineage/out.ipynb",
    params:
        dom_indir = lambda wildcards, input: dirname(input.dominant_cluster),
        script = join(ROOT_DIR, "src", "clone_cluster_count_embed", "aggregate_clone_meta.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note)
    shell:
        "papermill -p dom_indir {params.dom_indir} -p enrichment_f {input.clone_change_f} -p outdir {params.outdir} {params.script} {output.note}"

rule merge_enrich_and_lineage_and_input:
    input:
        clone_change_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/clone_change.csv",
        dominant_cluster = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/dominant_clone_clust/combinedConditions_dominant_cluster.csv",
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/SE.rds",
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/merge_enrich_and_lineage_and_input/out.ipynb",
    params:
        dom_indir = lambda wildcards, input: dirname(input.dominant_cluster),
        script = join(ROOT_DIR, "src", "clone_cluster_count_embed", "aggregate_clone_meta.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note)
    shell:
        "papermill -p se_f {input.se_f} -p dom_indir {params.dom_indir} -p enrichment_f {input.clone_change_f} -p outdir {params.outdir} {params.script} {output.note}"


use rule runDE_TF from annoMod as anno_runDE_TF with:
    input:
        integrated="{out_dir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds"
    output:
        note="{out_dir}/data/annotation/gff_{gff}/mergedSamples/DE_TF/DE_TF.ipynb"


use rule runGSEA from annoMod as anno_runGSEA with:
    input:
        DE_out_path = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/DE/DE.ipynb"
    output:
        note = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/DE/GSEA/clusters/GSEA.ipynb"
