report: "report/workflow.rst"

#from src.config import ROOT_DIR
#workdir: ROOT_DIR
wildcard_constraints:
    cellr='True|False'

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


# Merge the experiment-specific config with the pipeline parameters
params["results"] = res
#params["samples_meta"] = config["samples_meta"]
params["samples"] = samples
params["genome"] = config["genome"]
params["N_DONORS"] = config["N_DONORS"]
params["anno_res"] = anno_res
config["samples"] = samples
config["params"] = params
####

params_mt = params["mtpreproc"]
params_ft = params["filters"]
ft = params_ft["params"]
params_mult = params["multiplex"]
params_mgatk = params["mgatk"]
params_clones = params["clones"]
params_annclo = params["annotation_clones"]["params"]
params_clch = params["de_clones_change"]["params"]
params_clch["samples"] = samples

cellrbc = params_mt["params"]["cellrbc"]
num_reads_filter = params_mt["params"]["numreadsfilter"]
maxBP = params_mt["maxBP"]
ref_fa = params["genome_path"][config["genome"]]["ref_fa"]
gff = params["genome_path"][config["genome"]]["gff"]


ncellsthreshmgatk = params_mgatk["params"]["ncellsthresh"]
num_cells = params_mult["pseudo_multiplex"]["num_cells"]
is_prop = params_mult["pseudo_multiplex"]["is_proportional"]

nclonelist = params_clones['vireo']['params']['nclonelist']
config["params_clones"] = params_clones


params["markers_f"] = config["markers_f"]

# Add cell markers to annotation:
cluster_names_f = join(res, "results", "atac_clusters.txt")
if not os.path.exists(cluster_names_f):
    cluster_names_f = "None"
params["cluster_names_f"] = cluster_names_f
#workdir: config["work_dir"]


###################################
## Rules and workflow
###################################
###################################
rule all:
    input:
        expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/thr__{t}_rt__{rt}/annotation_clones/de_clone_btwnvars_RNA_af/mtPlots.ipynb",
            output=res,cellrbc=cellrbc,num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            kparam=config["params_clones"]["knn"]["params"]["resolution"],
            variants=[x for x in params_clones["variants"] if x != "simple"], gff=gff,
            t=params_annclo["t"], rt=params_annclo["rt"]),
        # dominant cluster and umap clones
        expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/clone_clust_embed/tsne_perp{perp}_donperp{donperp}/embed.ipynb",
            output=res,cellrbc=cellrbc,num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            kparam=config["params_clones"]["knn"]["params"]["resolution"],
            variants=[x for x in params_clones["variants"] if x != "simple"], gff=gff,
            perp=params["clone_clust_embed"]["params"]["perplexity"], donperp=params["clone_clust_embed"]["params"]["donor_perplexity"]),
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/{de_group}/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/{dnr_grp}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
        #     output=res,cellrbc=cellrbc,num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"],
        #     variants=[x for x in params_clones["variants"] if x != "simple"],
        #     logfc_threshold=params_clch["logfc_threshold"], p_thresh=params_clch["pthresh"],
        #     minPct=params_clch["min_pct"], gsea_pval=params_clch["gsea_pval"], #assay=params_annclo["assay"],
        #     gff=gff, prefilter=params_clch["prefilter"],
        #     stat_col=params_clch["stat_col"], padjmethod=params_clch["padjmethod"],
        #     filt = params_clch["filt"], shuffle=params_clch["shuffle"], use_p_adjust=params_clch["use_p_adjust"],
        #     pthresh=params_clch["pthresh"], min_cell=params_clch["min_cell"], min_cell_both=params_clch["min_cell_both"],
        #     de_group=["btwnChange", "btwnChange_inClust"] , dnr_grp=["allDonors", "sepDonors"]
        #     ),
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/allDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
        #     output=res,cellrbc=cellrbc,num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"],
        #     variants=[x for x in params_clones["variants"] if x != "simple"],
        #     logfc_threshold=params_annclo["logfc_threshold"], p_thresh=params_annclo["p_thresh"],
        #     minPct=params_annclo["min_pct"], gsea_pval=params_annclo["gsea_pval"], #assay=params_annclo["assay"],
        #     gff=gff, min_cells=params_annclo["min_cells"], prefilter=params_annclo["prefilter"],
        #     stat_col=params_annclo["stat_col"], padjmethod=params_annclo["padjmethod"],
        #     filt = params_clch["filt"], shuffle=params_clch["shuffle"], use_p_adjust=params_clch["use_p_adjust"],
        #     pthresh=params_clch["pthresh"], min_cell=params_clch["min_cell"], min_cell_both=params_clch["min_cell"],
        # ),
        expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/thr__{t}_rt__{rt}/annotation_clones/de_clone_btwnvars_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
            output=res,cellrbc=cellrbc,num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            kparam=config["params_clones"]["knn"]["params"]["resolution"],
            variants=[x for x in params_clones["variants"] if x != "simple"],
            logfc_threshold=params_annclo["logfc_threshold"], p_thresh=params_annclo["p_thresh"],
            btwnMinpct=params_annclo["btwnMinpct"], gsea_pval=params_annclo["gsea_pval"], #assay=params_annclo["assay"],
            gff=gff, min_cells=params_annclo["min_cells"], prefilter=params_annclo["prefilter"],
            stat_col=params_annclo["stat_col"], padjmethod=params_annclo["padjmethod"],
            t=params_annclo["t"], rt=params_annclo["rt"]),

        expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/clone_counts/minCellConds_{min_cells}/clone_sizes_norm.csv",
            output=res,cellrbc=cellrbc,num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            kparam=config["params_clones"]["knn"]["params"]["resolution"],
            variants=[x for x in params_clones["variants"] if x != "simple"],
            gff=gff, min_cells=params_annclo["min_cells"]),
        expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/result.csv",
            output=res,cellrbc=cellrbc,num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            kparam=config["params_clones"]["knn"]["params"]["resolution"],
            variants=[x for x in params_clones["variants"] if x != "simple"],
            logfc_threshold=params_annclo["logfc_threshold"], p_thresh=params_annclo["p_thresh"],
            btwnMinpct=params_annclo["btwnMinpct"], gsea_pval=params_annclo["gsea_pval"], #assay=params_annclo["assay"],
            gff=gff, hyperMinCl=params_annclo["hyperMinCl"],bothConds=params_annclo["bothConds"], pthresh=params_annclo["p_thresh"]),
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/de_clone_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
        #     output=res,cellrbc=cellrbc,num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"],
        #     variants=[x for x in params_clones["variants"] if x != "simple"],
        #     logfc_threshold=params_annclo["logfc_threshold"], p_thresh=params_annclo["p_thresh"],
        #     btwnMinpct=params_annclo["btwnMinpct"], gsea_pval=params_annclo["gsea_pval"], #assay=params_annclo["assay"],
        #     gff=gff, min_cells=params_annclo["min_cells"], prefilter=params_annclo["prefilter"],
        #     stat_col=params_annclo["stat_col"], padjmethod=params_annclo["padjmethod"]),
        #
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/de_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_p{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
        #     output=res,cellrbc=cellrbc,num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"],
        #     variants=[x for x in params_clones["variants"] if x != "simple"],
        #     logfc_threshold=params_annclo["logfc_threshold"], p_thresh=params_annclo["p_thresh"],
        #     btwnMinpct=params_annclo["btwnMinpct"], gsea_pval=params_annclo["gsea_pval"], #assay=params_annclo["assay"],
        #     gff=gff, min_cells=params_annclo["min_cells"], prefilter=params_annclo["prefilter"],
        #     stat_col=params_annclo["stat_col"], padjmethod=params_annclo["padjmethod"]),

        expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/markers/markers.ipynb",
            output=res,cellrbc=cellrbc,num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            kparam=config["params_clones"]["knn"]["params"]["resolution"],
            variants=[x for x in params_clones["variants"] if x != "simple"],
            gff=gff),
        expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/de_btwnclust_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}/de_p{p_thresh}_GSEA_pthresh_{gsea_pval}/summary.ipynb",
            output=res,cellrbc=cellrbc,num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            kparam=config["params_clones"]["knn"]["params"]["resolution"],
            variants=[x for x in params_clones["variants"] if x != "simple"],
            logfc_threshold=params_annclo["logfc_threshold"], p_thresh=params_annclo["p_thresh"],
            btwnMinpct=params_annclo["btwnMinpct"], gsea_pval=params_annclo["gsea_pval"], #assay=params_annclo["assay"],
            gff=gff),
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/de_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_p{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
        #     output=res,cellrbc=cellrbc,num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"],
        #     variants=[x for x in params_clones["variants"] if x != "simple"],
        #     logfc_threshold=params_annclo["logfc_threshold"], p_thresh=params_annclo["p_thresh"],
        #     btwnMinpct=params_annclo["btwnMinpct"], #assay=params_annclo["assay"],
        #     gff=gff, gsea_pval=params_annclo["gsea_pval"], prefilter=params_annclo["prefilter"],
        #     stat_col=params_annclo["stat_col"], padjmethod=params_annclo["padjmethod"]),
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/{f}",
        #     f = ["dotplot.allDonors.clones.topPvals.pdf","dotplot.allDonors.clones.topPvals.noZ.pdf", "clone.TFactivity.topPvals.txt"],
        #     output=res,cellrbc=cellrbc,num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"],
        #     variants=[x for x in params_clones["variants"] if x != "simple"],
        #     logThresh=params["annotation_clones"]["params"]["logfc_threshold"],
        #     minPct=params["annotation_clones"]["params"]["min_pct"],
        #     cdf_thresh=params["annotation_clones"]["params"]["cdf_thresh"]),
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/{f}",
        #     f=["dotplot.allDonors.clones.topPvals.pdf","dotplot.allDonors.clones.topPvals.noZ.pdf", "clone.TFactivity.topPvals.txt"],
        #     output=res,cellrbc=cellrbc,num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft[
        #         'minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
        #     nclones=nclonelist,
        #     variants=params_clones["variants"],
        #     logThresh=params["annotation_clones"]["params"]["logfc_threshold"],
        #     minPct=params["annotation_clones"]["params"]["min_pct"],
        #     cdf_thresh=params["annotation_clones"]["params"]["cdf_thresh"]),

        # A. Preproc: Initial MT coverage per position per cell
        expand("{output}/figures/{sample}/MT/cellr_{cellrbc}/numread_{num_read}_MT_position_coverage.png",  output=res,sample=samples.index, num_read=num_reads_filter, cellrbc=cellrbc),
         # B. Multiplex: Multiplexed VAF and depth dendrograms for each donor and selected variants
        # expand("{output}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/donor{d}.pdf",
        #        output=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #        mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #        hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], d=np.arange(config["N_DONORS"])),
        expand("{output}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/multiplex.pdf",
               output=res,cellrbc=cellrbc, num_read=num_reads_filter,
               mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
               hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh']),

        ## Clones
        expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/{method}/.completed",
            output=res,cellrbc=cellrbc,num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],
            minhetcells=ft['minhetcells'], hetcountthresh=ft['hetcountthresh'],
            bqthresh=ft['bqthresh'],method=params_clones["method"]),
        #
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/.pipeline",
        #     output=res, cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], method=params_clones["method"], variants=params_clones["variants"]),
        # #vireo
        # expand("{output}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/clones_dendro/donor{d}.na.clust.max2.AF.png",
        #     output=res,cellrbc=cellrbc,num_read=num_reads_filter,
        #     d=np.arange(config['N_DONORS']),
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
        #     method=params_clones["method"], nclones=nclonelist,
        #     variants=params_clones["variants"]
        # ),
        # # knn
        # expand(
        #     "{output}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/clones_dendro/donor{d}.na.clust.max2.AF.png",
        #     output=res, cellrbc=cellrbc, num_read=num_reads_filter,
        #     d=np.arange(config['N_DONORS']),
        #         mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft[
        #             'minhetcells'],
        #         hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
        #         method=params_clones["method"],
        #         kparam=params_clones["knn"]["params"]["resolution"],
        #         variants=[x for x in params_clones["variants"] if x!="simple"]
        #     ),
        expand(
            "{output}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/clones_barcodes/donor{d}.ipynb",
            output=res, cellrbc=cellrbc, num_read=num_reads_filter,
            d=np.arange(config['N_DONORS']),
                mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft[
                    'minhetcells'],
                hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
                method=params_clones["method"],
                kparam=params_clones["knn"]["params"]["resolution"],
                variants=[x for x in params_clones["variants"] if x!="simple"]
            ),
        ## Enrichment shuffle stats
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/enrichment/shuffle_stats/shuffle_stats.csv",
        #     output=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
        #     variants=params_clones["variants"], nclones=nclonelist),
        expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/shuffle_stats/shuffle_stats.csv",
            output=res,cellrbc=cellrbc, num_read=num_reads_filter,
            variants=[x for x in params_clones["variants"] if x != "simple"],
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
            kparam=params_clones["knn"]["params"]["resolution"]),
        ##  Methods compare
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/comparisons/comparisons.ipynb",
        #     output=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh']),

        # # C. Clones: Variant types
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb",
        #        output=res, cellrbc=cellrbc, num_read=num_reads_filter,
        #        mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #        hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], n_clones=config['multiplex']['n_clone_list']),

        # # D. Enrichment: Volcano plot
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/enrichment/.status",
        #        output=res, cellrbc=cellrbc, num_read=num_reads_filter,
        #        mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #        hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh']),
        # # E. Lineage labels for cells, which can be input to loupe browser
        # expand("{output}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/clones/variants_{variants}/{method}_{params_name}{method_params}/donor{d}/cells_BC.csv",
        #        d=range(config["N_DONORS"]),
        #        method=config["clones"]["method"],
        #        method_params=[config[x]['params'] for x in config["method"]])
        #

################################################################
## Import from prior snakefile modules
## Here, we redefine the input to be based on our config['files'] dictionary
from snakemake.utils import min_version
min_version("6.0")
