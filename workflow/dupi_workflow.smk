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
#anno_res = join(config["outdir"], "annotation", "data", config["prefix"])
anno_res = join(res, "data", "annotation")


# Merge the experiment-specific config with the pipeline parameters
params["results"] = res
#params["samples_meta"] = config["samples_meta"]
params["samples"] = samples
params["genome"] = config["genome"]
params["N_DONORS"] = config["N_DONORS"]
params["anno_res"] = anno_res
config["samples"] = samples
config["params"] = params
config["results"] = res
####

params_mt = params["mtpreproc"]
params_ft = params["filters"]
ft = params_ft["params"]
params_mult = params["multiplex"]
params_mgatk = params["mgatk"]
params_clones = params["clones"] # Not params b/c there are separate methods here
params_annclo = params["annotation_clones"]["params"]
params_clch = params["de_clones_change"]["params"]
params_clch["samples"] = samples

cellrbc = params_mt["params"]["cellrbc"]
num_reads_filter = params_mt["params"]["numreadsfilter"]
maxBP = params_mt["maxBP"]
ref_fa = params["genome_path"][config["genome"]]["ref_fa"]
mt_ref_fa = params["genome_path"][config["genome"]]["mt_ref_fa"]
config["mt_ref"] = mt_ref_fa

gff = params["genome_path"][config["genome"]]["gff"]

params["gff"] = gff

ncellsthreshmgatk = params_mgatk["params"]["ncellsthresh"]
num_cells = params_mult["pseudo_multiplex"]["num_cells"]
is_prop = params_mult["pseudo_multiplex"]["is_proportional"]

nclonelist = params_clones['vireo']['params']['nclonelist']
config["params_clones"] = params_clones


params["markers_f"] = config["markers_f"]

# Add cell markers to annotation:
cluster_names_f = config.get("umap_clusters_f", None) #join(res, "results", "atac_clusters.txt")
if (cluster_names_f is not None) and (not os.path.exists(cluster_names_f)):
    print("Cluster label file not here.")
    cluster_names_f = ""
params["cluster_names_f"] = cluster_names_f
#workdir: config["work_dir"]


# somvars parameters
if "somvars_dir" in config:
    somvar_res = join(config["outdir"], "somatic_variants", "data", config["somvars_dir"])
else:
    somvar_res = join(config["outdir"], "somatic_variants", "data", "somvars")

peak_names = [] if "bed_regions" not in config else list(config["bed_regions"].keys())
peak_names.append("all")
print("peak_names")

gatk_peak_names = [] if "gatk_bed_regions" not in config else list(config["gatk_bed_regions"].keys())

if "som_vars" in params.keys():
    print(params["som_vars"]["mutect"])
    params_somvar = params["som_vars"]["mutect"]["params"]
else:
    params_somvar = {}

###################################
## Rules and workflow
###################################
###################################
rule all:
    input:
        # # A. Preproc: Initial MT coverage per position per cell
        # expand("{out_dir}/figures/{sample}/MT/cellr_{cellrbc}/numread_{num_read}_MT_position_coverage.png",
        #     out_dir=res,sample=samples.index, num_read=num_reads_filter, cellrbc=cellrbc),
        # B. Multiplex: Multiplexed VAF and depth dendrograms for each donor and selected variants
        expand("{out_dir}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/multiplex.pdf",
            out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh']),
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/elbo/out_multiplex_elbo.ipynb",
            out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            method=params_clones["method"],
            kparam=params_clones["knn"]["params"]["resolution"],d=np.arange(config["N_DONORS"])),

         #multiplex variant types
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_{variants}/variantTypes/output.ipynb",
            out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            variants=[x for x in params_clones["variants"] if x!="simple"]),

        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_init/donor{d}/af.tsv",
            out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            method=params_clones["method"],
            kparam=params_clones["knn"]["params"]["resolution"],d=np.arange(config["N_DONORS"])),

        # Variants
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_variants/anno_variants.tsv",
            out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh']),

        # Clones (lineage)
        # donor specific
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
        #     out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'], variants=[x for x in params_clones["variants"] if x!="simple"],
        #     kparam=params_clones["knn"]["params"]["resolution"],d=np.arange(config["N_DONORS"])),
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_init/knn/kparam_{kparam}/cells_meta.tsv",
        #     out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
        #     method=params_clones["method"],
        #     kparam=params_clones["knn"]["params"]["resolution"],
            #variants=[x for x in params_clones["variants"] if x!="simple"]
       # ),
        # Clone stats
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/_clone_complete.txt",
            out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            method=params_clones["method"],
            kparam=params_clones["knn"]["params"]["resolution"],
            variants=[x for x in params_clones["variants"] if x!="simple"]
        ),

        # # Enrichment
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/_enrichment_complete",
            out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
            variants=[x for x in params_clones["variants"] if x != "simple"],
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
            kparam=params_clones["knn"]["params"]["resolution"]),

        # Annotation clones
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/_nuclear_clones_complete.txt",
            out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], gff=gff,
            variants=[x for x in params_clones["variants"] if x != "simple"],
            kparam=params_clones["knn"]["params"]["resolution"]),

        # MT Plots
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/thr__{t}_rt__{rt}/annotation_clones/de_clone_btwnvars_RNA_af/mtPlots.ipynb",
        #     out_dir=res,cellrbc=cellrbc,num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"],
        #     variants=[x for x in params_clones["variants"] if x != "simple"], gff=gff,
        #     t=params_annclo["t"], rt=params_annclo["rt"]),

        # Cluster DE
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/de_btwnclust_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}/de_p{p_thresh}_GSEA_pthresh_{gsea_pval}/summary.ipynb",
            out_dir=res,cellrbc=cellrbc,num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            kparam=config["params_clones"]["knn"]["params"]["resolution"],
            variants=[x for x in params_clones["variants"] if x != "simple"], gff=gff,
            p_thresh=params_annclo["p_thresh"],
            min_cells=params_annclo["min_cells"],
            logfc_threshold=params_annclo["logfc_threshold"],
            btwnMinpct=params_annclo["min_pct"], gsea_pval=params_annclo["gsea_pval"],
            ),


        #
        # # Condition DE
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/de_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_p{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
            out_dir=res,cellrbc=cellrbc,num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            kparam=config["params_clones"]["knn"]["params"]["resolution"],
            variants=[x for x in params_clones["variants"] if x != "simple"], gff=gff,
            p_thresh=params_annclo["p_thresh"],
            logfc_threshold=params_annclo["logfc_threshold"],
            btwnMinpct=params_annclo["min_pct"], gsea_pval=params_annclo["gsea_pval"],
            prefilter=params_annclo["prefilter"],
            stat_col=params_annclo["stat_col"], padjmethod=params_annclo["padjmethod"]),

        #
        # # Within each clone run DE btwn conditions
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/de_clone_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
        #     out_dir=res,cellrbc=cellrbc,num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"],
        #     variants=[x for x in params_clones["variants"] if x != "simple"], gff=gff,
        #     p_thresh=params_annclo["p_thresh"],
        #     min_cells=params_annclo["min_cells"],
        #     logfc_threshold=params_annclo["logfc_threshold"],
        #     btwnMinpct=params_annclo["min_pct"], gsea_pval=params_annclo["gsea_pval"],
        #     prefilter=params_annclo["prefilter"],
        #     stat_col=params_annclo["stat_col"], padjmethod=params_annclo["padjmethod"]),
        #
        # # Clone expand-regress DE
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/clones_change/_de_clone_change_complete.txt",
        #             out_dir=res,cellrbc=cellrbc,num_read=num_reads_filter,
        #             mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #             hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
        #             kparam=config["params_clones"]["knn"]["params"]["resolution"],
        #             variants=[x for x in params_clones["variants"] if x != "simple"], gff=gff,
        #             ),

        ##  Methods compare
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/comparisons/comparisons.ipynb",
            out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh']),

        # # A. Merge clone enrich and nuclear B. mt as clones
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/merge_enrich_and_lineage/out.ipynb",
        #         out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #         mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #         hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
        #         kparam=config["params_clones"]["knn"]["params"]["resolution"],
        #         minPct=params_annclo["min_pct"], gsea_pval=params_annclo["gsea_pval"], #assay=params_annclo["assay"],
        #         gff=gff, min_cells=params_annclo["min_cells"],
        #         variants=[x for x in params_clones["variants"] if x != "simple"],
        #         stat_col=params_annclo["stat_col"], padjmethod=params_annclo["padjmethod"],
        #         filt = params_clch["filt"], shuffle=params_clch["shuffle"], use_p_adjust=params_clch["use_p_adjust"],
        #         pthresh=params_clch["pthresh"], min_cell=params_clch["min_cell"], min_cell_both=params_clch["min_cell_both"],),

        # mt as clones
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/mt_as_clones/variants_{variants}/.complete.txt",
                out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
                mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
                hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], variants=[x for x in params_clones["variants"] if x != "simple"],),
        #
        # # clonal_shifts
        #  expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clonal_shifts/variants_{variants}/.finalize_complete.txt",
        #         out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #         mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #         hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], variants=[x for x in params_clones["variants"] if x != "simple"]),

        # # individual clones
        #  expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/single_clones/{f}",
        #         out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #         mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #         hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], f = [".preproc.txt", ".top_clones.txt"]) #".aggregate.txt"

        #"{outdir}/single_clones/.preproc.txt"

        # # C. Clones: Variant types



        # Somatic variants

        #clones_variants_summary/summary_som.ipynb
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/gatk_mutect/peaks_{regions}/minaf{minaf}_midaf{midaf}_ad{minad}_dp{mindp}_qual{qual}/som_dendro_{d_thresh}/clones_variants_summary/summary_som.ipynb",
        #     out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"], variants=[x for x in params_clones["variants"] if x!="simple"], d_thresh=params_clones["params"]["dendro_thresh"],
        #     sample=samples.index, regions=gatk_peak_names,
        #     minaf=params_somvar["min_af"], midaf=params_somvar["mid_af"], minad=params_somvar["min_ad"],
        #     mindp=params_somvar["min_dp"], qual=params_somvar["qual"]),
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/gatk_mutect/peaks_{regions}/minaf{minaf}_midaf{midaf}_ad{minad}_dp{mindp}_qual{qual}/som_dendro_{d_thresh}/som_clones.ipynb",
        #     out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"], variants=[x for x in params_clones["variants"] if x!="simple"], d_thresh=params_clones["params"]["dendro_thresh"],
        #     sample=samples.index, regions=gatk_peak_names,
        #     minaf=params_somvar["min_af"], midaf=params_somvar["mid_af"], minad=params_somvar["min_ad"],
        #     mindp=params_somvar["min_dp"], qual=params_somvar["qual"]),
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/gatk_mutect/peaks_{regions}/minaf{minaf}_midaf{midaf}_ad{minad}_dp{mindp}_qual{qual}/som_dendro_{d_thresh}/filt_chip/som_clones.ipynb",
        #     out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"], variants=[x for x in params_clones["variants"] if x!="simple"], d_thresh=params_clones["params"]["dendro_thresh"],
        #     sample=samples.index, regions=gatk_peak_names,
        #     minaf=params_somvar["min_af"], midaf=params_somvar["mid_af"], minad=params_somvar["min_ad"],
        #     mindp=params_somvar["min_dp"], qual=params_somvar["qual"]),
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/needle/peaks_{regions}/{sample}/som_clones.ipynb",
        #     out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"], variants=[x for x in params_clones["variants"] if x!="simple"],
        #     sample=samples.index, regions=peak_names),
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/needle/peaks_{regions}/som_dendro_{d_thresh}/som_clones.ipynb",
        #     out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"], variants=[x for x in params_clones["variants"] if x!="simple"], d_thresh=params_clones["params"]["dendro_thresh"],
        #     sample=samples.index, regions=peak_names),
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/needle/peaks_{regions}/som_dendro_{d_thresh}/clones_variants_summary/summary_som.ipynb",
        #     out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"], variants=[x for x in params_clones["variants"] if x!="simple"], d_thresh=params_clones["params"]["dendro_thresh"],
        #     sample=samples.index, regions=peak_names,
        #     minaf=params_somvar["min_af"], midaf=params_somvar["mid_af"], minad=params_somvar["min_ad"],
        #     mindp=params_somvar["min_dp"], qual=params_somvar["qual"]),

        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/gatk_mutect/peaks_{regions}/som_dendro_{d_thresh}/chip_clones_variants_summary/summary_som.ipynb",
        #     out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
        #     kparam=config["params_clones"]["knn"]["params"]["resolution"], variants=[x for x in params_clones["variants"] if x!="simple"], d_thresh=params_clones["params"]["dendro_thresh"],
        #     sample=samples.index, regions=gatk_peak_names),

#"{outdir}/somatic_variants/{somvar_method}/peaks_{regions}/filt_chip/som_dendro_{dendro_thresh}/som_clones.ipynb",
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb",
        #        output=res, cellrbc=cellrbc, num_read=num_reads_filter,
        #        mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #        hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], n_clones=config['multiplex']['n_clone_list']),

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




########################################################################
########################################################################


########################################################################
## 1. Go from 10x output to filtered scPileup outdir.
########################################################################
use rule * from mtpreprocMod as mtpreproc_*

########################################################################
## 2. Call variants using MGATK and convert to the Vireo multiplex input
########################################################################
use rule * from mgatkMod

use rule mgatk_to_vireoIn from mgatkMod with:
    input:
        "{out_dir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/{sample}.variant.rds"
    output:
        "{out_dir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/{sample}/vireoIn/cellSNP.tag.AD.mtx"


rule mgatk_to_variant:
    input:
        cells = lambda wildcards: expand("{{out_dir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/{sample}/vireoIn/cellSNP.tag.AD.mtx",
                                  sample=samples.index)
    output:
        vcf = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_variants/variants.vcf",
    params:
        note = join(ROOT_DIR, "workflow", "notebooks", "variant_types", "scoreVariants.ipynb"),
        vcfs = lambda wildcards, input: [join(dirname(x),"cellSNP.base.vcf") for x in input.cells],
        mt_fasta = mt_ref_fa
    run:
        #all_vcf = pd.read_csv(params.vcfs[0], sep="\t")[["#CHROM", "POS", "REF", "ALT"]]
        #for i in params.vcfs[1:]:
            #all_vcf = all_vcf.merge(pd.read_csv(i, sep="\t")[["#CHROM", "POS", "REF", "ALT"]], how="inner")
        all_vcf = []
        for i in params.vcfs:
            all_vcf.append(pd.read_csv(i, sep="\t")[["#CHROM", "POS", "REF", "ALT"]])
        all_vcf = pd.concat(all_vcf, ignore_index=True).drop_duplicates()
        all_vcf["QUAL"] = "."
        all_vcf["FILTER"] = "."
        all_vcf["INFO"] = "."
        all_vcf["ID"] =  all_vcf["REF"] +  ">" + all_vcf["ALT"]
        all_vcf["REF"] = [x[-1] for x in all_vcf["REF"]]
        header = "##fileformat=VCFv4.0"
        header = header + "\n" + f"##reference=file:/{params.mt_fasta}\n"
        print('header', header)
        print('output.vcf', output.vcf)
        with open(output.vcf, 'a') as file:
            file.write(header)
        print(all_vcf.columns)
        all_vcf.to_csv(output.vcf, mode='a', index=False, sep="\t", header=True)


rule merge_variant_samples:
    input:
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_variants/variants.vcf",
    output:
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_variants/anno_variants.tsv",
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_variants/anno_variants.ipynb"
    params:
        note = join(ROOT_DIR, "workflow", "notebooks", "variant_types", "scoreVariants.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
    shell:
        "papermill -p vcf {input[0]} -p outdir {params.outdir} {params.note} {output.note}"

########################################################################
## 3. Merge the sample inputs together by concatenating.
########################################################################
##### The cell barcodes get changed here!
rule merged:
    """ Merge the sample pileup matrices
    """
    input:
        lambda wildcards: expand("{{out_dir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/{sample}/vireoIn/cellSNP.tag.AD.mtx",
                                        sample=samples.index)# , minC=wildcards.minC, minAF=wildcards.minAF)
    output: "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/cellSNP.tag.AD.mtx"
    #input: lambda wildcards: expand("data/{{prefix}}/chrM/{name}_cellSNP_minC{{mt_minC}}_minAF{{mt_minAF}}", name=config["samples"])
    params:
        num_cells=num_cells,
        is_prop = is_prop,
        indirs = lambda wildcards, input: [dirname(x) for x in input],
        outdir = lambda wildcards, output: dirname(output[0]),
        prefix= ','.join(samples.index)
    #log: "logs/{prefix}/vireo/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}.log"
    shell: "python -m src.pseudo_batch {params.outdir} {params.indirs} --num_cells {params.num_cells} --is_prop {params.is_prop} --samples {params.prefix}" # > {log} 2>&1"


########################################################################
## 4. De-multiplex the merged conditions using Vireo
########################################################################
use rule * from multMod as mult_*

use rule donors_type_variants from multMod with:
    output:
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/variants/variants.ipynb",


use rule donors_plotAF from multMod as mult_donors_plotAF with:
    input:
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb",
    output:
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/dendrograms/figures/donor{d}_dendrogram.{suf}",
               d=np.arange(config["N_DONORS"]), suf=["png", "depth.png", "withHigh.png"]),
        #"{out_dir}/merged/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/dendrograms/af_dendro.ipynb",

rule multiplex_report:
    input:
        multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/",
                 "multiplex_AF_SNPs_all_afFilt.png", "multiplex_clusters_all.labels.png")
    output:
        report("{out_dir}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/multiplex.pdf", category="Multiplex")
    shell:
        "convert {input} {output[0]}"

rule donor_report:
    input:
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/dendrograms/figures/donor{{d}}_dendrogram.{suf}",
               suf=["png", "depth.png", "withHigh.png"]),
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/mgatk/d{d}.variantQC.png",
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/donor{{d}}.labels.png",
               nclones=nclonelist, variants=params_clones["variants"]),
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/donor{{d}}.variants.labels.png",
               nclones=nclonelist, variants=params_clones["variants"])
    output:
        report("{out_dir}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/donor{d}.pdf", category="Donor-specific")
    shell:
        "convert {input} {output[0]}"


rule donors_se_meta:
    """Save seurat meta-data as csv"""
    input:
        se_f = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds",
        cells_meta_dir = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex_AF_SNPs_all_afFilt.png"
    output:
        se_meta = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_multiplex/gff_{gff}/se_cells_meta.tsv",
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_multiplex/gff_{gff}/se_cells_meta.ipynb",
        se = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_multiplex/gff_{gff}/SE.rds",
    params:
        outdir = lambda wildcards, output: dirname(output.note),
        cells_meta = lambda wildcards, input: join(dirname(input.cells_meta_dir), "cells_meta.tsv"),
        rscript = join(ROOT_DIR, "workflow/notebooks/lineage_clones/addDonors.ipynb"),
    shell: "papermill -p se_f {input.se_f} -p outdir {params.outdir} -p cells_meta_f {params.cells_meta}  -p save_rds TRUE {params.rscript} {output.note}"

def get_cluster_labels():
    out = config.get("umap_clusters_f", "FALSE")
    if out is None :
        out = "FALSE"
    return out

rule don_add_cluster_labels:
    """Prepare clone-by-cluster counts for umap and hypergeometric test"""
    input:
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_multiplex/gff_{gff}/se_cells_meta.tsv",
    output:
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_multiplex/gff_{gff}/se_cells_meta_labels.tsv",
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_multiplex/gff_{gff}/add_cluster_labels.ipynb",
    params:
        labels = get_cluster_labels(),
        rscript = join(ROOT_DIR, "workflow/notebooks/lineage_clones/add_cluster_labels.ipynb"),
    shell: "papermill -p se_f {input.se_f} -p cluster_labels_f {params.labels} -p is_rds FALSE {params.rscript} {output.note}"


rule mt_to_clones:
    input:
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/clones_init/donor{d}/af.tsv",
                d=np.arange(config["N_DONORS"]))
    output:
        note="{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/mt_clones/thr__{t}_rt__{rt}/mtClones.ipynb",
        mtClones="{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/mt_clones/thr__{t}_rt__{rt}/cells_meta.tsv"
    params:
        don_dir = lambda wildcards, input: dirname(dirname(input[0])),
        outdir = lambda wildcards, output: dirname(output.note),
        rscript = join(ROOT_DIR, "src", "clones", "mt_to_clones.ipynb"),
        n_donors = config["N_DONORS"],
        samples = ",".join(samples['sample_name'].values)
    shell: "papermill -p don_dir {params.don_dir} -p outdir {params.outdir} -p t {wildcards.t} -p rt {wildcards.rt} -p n_donors {params.n_donors} -p samples {params.samples} {params.rscript} {output.note}"


########################################################################
## 5. Calling clones for each donor
########################################################################
rule vireo:
    input:
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb" #get_input_multiplex
    output:
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{{nclones}}/donor{d}.labels.png", d=np.arange(config["N_DONORS"])),#, category="lineage"),
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{{nclones}}/donor{d}.variants.labels.png", d=np.arange(config["N_DONORS"])),#, category="lineage"),
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{nclones}/cells_meta.tsv"
    params:
        notebook=join("src", "vireo", "2_MT_Lineage_Construct.ipynb"),
        #output_notebook = lambda wildcards, output: join(dirname(output[0]), 'clones.ipynb'),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["N_DONORS"],
        nclones= lambda wildcards: wildcards.nclones #",".join([str(x) for x in nclonelist])
    threads: 8
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p nclones {params.nclones} {params.notebook} {params.OUTDIR}/output.ipynb"


# ########################################################################
# ## Variants A: After multiplexing, separate by donor, grab variants from filters
# ## that overlap with both conditions, and then run mgatk to call variants again.
# ##
# # Extract pileups for each donor from original filter matrices
# ########################################################################
rule scPileup_filter_mgatk:
    input:
        cov =  expand("{{out_dir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/{sample}.coverage.txt",
            sample=samples.index),
        mult = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb",
    output:
        cov = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/d{d}.coverage.txt",
            d=range(config["N_DONORS"]))
    params:
        cells_meta = lambda wildcards, input: join(dirname(input.mult), "cells_meta.tsv"),
        sample = samples['sample_name'].values,
        #outdir = lambda wildcards, output: dirname(output.cov),
    script: join(ROOT_DIR, "src/donor_filter_mgatk.py")


module procDonMgatkMod:
    snakefile: "./rules/clone_detection/proc_donormgatk.smk"
    config: config

use rule * from procDonMgatkMod as donMGATK_*

# use rule knn_mgatkdonor_concat from procDonMgatkMod as donMGATK_knn_mgatkdonor_concat with:
#     input:
#         cells = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{{kparam}}/labels/donor{d}/cells_meta.tsv",
#                 d = np.arange(config["N_DONORS"])),
#     output:
#         cells_all = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{kparam}/cells_meta.tsv",
# #         expand("{{outdir}}/clones/variants_mgatkdonor/vireo/nclones{{nclones}}/donor{d}_cells_meta.tsv",
# #                d = np.arange(config["N_DONORS"])),

# ########################################################################
# ## Variants B: After multiplexing, separate by donor, keep variants from filters
# ## that overlap with both conditions, and then use those for clonal calling (may include a lot of variants).
# ##
# # Extract pileups for each donor from original filter matrices
# ########################################################################
# rule scPileup_init:
#     input:
#         cov =  expand("{{out_dir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/{sample}.coverage.txt",
#             sample=samples.index),
#         mult = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb"
#     output:
#         af = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/clones_init/donor{d}/af.tsv",
#                     d=np.arange(config["N_DONORS"]))
#     params:
#         sample = samples['sample_name'].values,
#         cells_meta = lambda wildcards, input: join(dirname(input.mult),"cells_meta.tsv"),
#         ref_mt = params_mgatk["chrM_refAllele"]
#     script: join(ROOT_DIR, "src/donor_to_clonesInput/donor_to_af.py")


rule scPileup_init_v02:
    """ The simpleUnion version was not working properly. This script works for both of them"""
    input:
        cov =  expand("{{out_dir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/{sample}.coverage.txt",
            sample=samples.index),
        mult = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb"
    output:
        af = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/clones_init/donor{d}/af.tsv",
                    d=np.arange(config["N_DONORS"]))
    params:
        sample = samples['sample_name'].values,
        cells_meta = lambda wildcards, input: join(dirname(input.mult),"cells_meta.tsv"),
        ref_mt = params_mgatk["chrM_refAllele"],
        var_type = "init"
    version: "v02"
    notebook:
        join(ROOT_DIR, "workflow/notebooks/af_filter/donor_to_af_v02.ipynb")


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
# rule scPileup_simpleUnion:
#     input:
#         cov =  expand("{{out_dir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/{sample}.coverage.txt",
#             sample=samples.index),
#         mult = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb"
#     output:
#         af = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/clones_simpleUnion/donor{d}/af.tsv",
#             d=np.arange(config["N_DONORS"]))
#     params:
#         sample = samples['sample_name'].values,
#         cells_meta = lambda wildcards, input: join(dirname(input.mult),"cells_meta.tsv"),
#         ref_mt = params_mgatk["chrM_refAllele"]
#     script: join(ROOT_DIR, "src/donor_to_clonesInput/donor_to_af_simpleUnion.py")

rule scPileup_simpleUnion_v02:
    input:
        cov =  expand("{{out_dir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/{sample}.coverage.txt",
            sample=samples.index),
        mult = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb"
    output:
        af = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/clones_simpleUnion/donor{d}/af.tsv",
            d=np.arange(config["N_DONORS"]))
    params:
        sample = samples['sample_name'].values,
        cells_meta = lambda wildcards, input: join(dirname(input.mult),"cells_meta.tsv"),
        ref_mt = params_mgatk["chrM_refAllele"],
        var_type = "simpleUnion"
    version: "v02"
    notebook:
        join(ROOT_DIR, "workflow/notebooks/af_filter/donor_to_af_v02.ipynb")


rule scPileup_variants:
    input:
        af = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/clones_{{variants}}/donor{d}/af.tsv",
            d=np.arange(config["N_DONORS"]))
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_{variants}/variantTypes/output.ipynb"
    params:
        indir = lambda wildcards, input: dirname(dirname(input.af[0])),
        outdir = lambda wildcards, output: dirname(output.note),
        note = join(ROOT_DIR, "workflow/notebooks/variant_types/Donors_variantTypes.ipynb")
    shell: "papermill -p indir {params.indir} -p outdir {params.outdir} {params.note} {output.note}"

module procSimUnMod:
    snakefile: "./rules/clone_detection/proc_simpleUnion.smk"
    config: config

use rule * from procSimUnMod as procSimUnMod_*

#########################################################################
## Clones stats and QC
#########################################################################
def get_counts_in(wildcards):
    w = wildcards
    print('w', w)
    if w.variants == "mgatkdonor":
        return f"{w.out_dir}/data/merged/MT/cellr_{w.cellrbc}/numread_{w.num_read}/filters/minC{w.mincells}_minR{w.minreads}_topN{w.topN}_hetT{w.hetthresh}_hetC{w.minhetcells}_hetCount{w.hetcountthresh}_bq{w.bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{w.d}/mgatk/d{w.d}.variant.rds"
    elif w.variants == "simple":
        return f"{w.out_dir}/data/merged/MT/cellr_{w.cellrbc}/numread_{w.num_read}/filters/minC{w.mincells}_minR{w.minreads}_topN{w.topN}_hetT{w.hetthresh}_hetC{w.minhetcells}_hetCount{w.hetcountthresh}_bq{w.bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb"
    elif w.variants == "init":
        return f"{w.out_dir}/data/merged/MT/cellr_{w.cellrbc}/numread_{w.num_read}/filters/minC{w.mincells}_minR{w.minreads}_topN{w.topN}_hetT{w.hetthresh}_hetC{w.minhetcells}_hetCount{w.hetcountthresh}_bq{w.bqthresh}/mgatk/vireoIn/multiplex/clones_init/donor{w.d}/af.tsv"
    elif w.variants == "simpleUnion":
        return f"{w.out_dir}/data/merged/MT/cellr_{w.cellrbc}/numread_{w.num_read}/filters/minC{w.mincells}_minR{w.minreads}_topN{w.topN}_hetT{w.hetthresh}_hetC{w.minhetcells}_hetCount{w.hetcountthresh}_bq{w.bqthresh}/mgatk/vireoIn/multiplex/clones_simpleUnion/donor{w.d}/af.tsv"



if "dendro_thresh" not in config:
    config["dendro_thresh"] = params_clones["params"]["dendro_thresh"]


module cloneMod:
    snakefile: "./rules/clone_detection/clone.smk"
    config: config


use rule * from cloneMod as clone_*


use rule convert_to_af from cloneMod as clone_convert_to_af with:
    input:
        counts_in = get_counts_in,
        cells_meta = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/sc_af/donor{d}/sc_af.ipynb",
    params:
        var_type = lambda wildcards: wildcards.variants,
        note = join(ROOT_DIR, "workflow", "notebooks", "clone_af_dendrograms", "Convert_to_AF.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
        indir = lambda wildcards, input: dirname(input.cells_meta),
        counts_indir = lambda wildcards, input: dirname(input.counts_in),


# use rule barcodes_inClones from cloneMod as clone_barcodes_inClones with:
#     input:
#         counts_in = get_counts_in,
#         cells_meta = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/btwnClones/donor{d}.ipynb",
#     params:
#         var_type = lambda wildcards: wildcards.variants,
#         note=join(ROOT_DIR, "worklow/notebooks/clone_af_dendrograms", "MT_barcodes/inClones.ipynb"),
#         outdir=lambda wildcards, output: dirname(output.note),
#         indir= lambda wildcards, input: dirname(input.cells_meta),
#         counts_indir = lambda wildcards, input: dirname(input.counts_in)


########################################################################
## Enrichment
########################################################################
module cloneEnrichMod:
    snakefile: "./rules/clone_detection/enrichment.smk"
    config: params

use rule * from cloneEnrichMod as cloneEnrich_*

# rule vireo_enrichment:
#     version: '1.0' # no +1 norm error
#     input:
#         "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/cells_meta.tsv"
#     output:
#         report("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/enrichment/volcano_Fisher_foldNorm.png",
#             category="lineage", subcategory="enrichment"),
#         report("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/enrichment/enrichmentNorm.csv",
#             category="lineage", subcategory="enrichment")
#     params:
#         #clones_indir = lambda wildcards, input: dirname(input[0]),#lambda wildcards, input: dirname(input[0]),
#         OUTDIR = lambda wildcards, output: dirname(output[0]),
#         script = join("src", "lineage", "lineage_enrichment.py"),
#         samples=",".join(samples.index),
#         comparisons = config["comparisons"] if "comparisons" in config else "None"
#     shell: "python {params.script} {input} {params.OUTDIR} {params.samples} --tests {params.comparisons:q}"
#
#
#
# rule knn_enrichment:
#     input:
#         "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
#     output:
#         report("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/volcano_Fisher_foldNorm.png",
#             category="lineage", subcategory="enrichment"),
#         report("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/enrichmentNorm.csv",
#             category="lineage", subcategory="enrichment")
#     params:
#         OUTDIR = lambda wildcards, output: dirname(output[0]),
#         script = join("src", "lineage", "lineage_enrichment.py"),
#         samples = ",".join(samples.index),
#         comparisons = config["comparisons"] if "comparisons" in config else "None"
#     shell: "python {params.script} {input} {params.OUTDIR} {params.samples} --tests {params.comparisons:q}"



print('params_clones', params_clones)
############################################
## Compare methods:
############################################
rule compare_methods:
    input:
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
                variants=[x for x in params_clones["variants"] if x != "simple"], kparam=params_clones["knn"]["params"]["resolution"]),
        # expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/cells_meta.tsv",
        #         variants=params_clones["variants"], nclones=nclonelist)
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/comparisons/comparisons.ipynb",
        fig = report(multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/comparisons/",
                              "methods_nAgreeNorm.png", "methods_nAgreeNorm_agg.png", "methods_n_T_T_agg.png", "methods_nTTNorm_agg.png"), category="Methods"),
    params:
        all = lambda wildcards, input: ",".join(x for x in input),
        note = join("src", "clones_compare", "distance_matrix.ipynb"),
        outdir = lambda wildcards, output: dirname(output[0]),
        prefix = config["prefix"],
        cfg_outdir = config["outdir"],
        params_f = config["config"]
    shell:
        "papermill -p all_files {params.all} -p outdir {params.outdir} -p cfg_outdir {params.cfg_outdir} -p prefix {params.prefix} -p params_f {params.params_f} {params.note} {output.note}"



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
params["umap_clusters_f"] = config.get("umap_clusters_f", "FALSE")

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


############################################
## clone expansion vs regressed DE
############################################
module cloneChangeMod:
    snakefile: "./rules/de/de_clones_change.smk"
    config: params #params_clch

use rule * from cloneChangeMod as cloneChange_*

use rule preprocess from cloneChangeMod as cloneChange_preprocess with:
    input:
        enrichment_dir = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/_enrichment_complete"
        #enrichment_dir = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/enrichmentNorm.csv"
    output:
        clone_change_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/clone_change.csv",
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/preproc.ipynb"

# use rule btwnClone_DE from cloneChangeMod as cloneChange_btwnClone_DE with:
#     input:
#         se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/SE.rds",
#         clone_change_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/clone_change.csv",
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/{de_group}/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{pthresh}/de.ipynb",


## get clone meta from stats and immune clusters
rule merge_enrich_and_lineage:
    input:
        clone_change_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/clone_change.csv",
        dominant_cluster = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/dominant_clone_clust/combinedConditions_dominant_cluster.csv"
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/merge_enrich_and_lineage/out.ipynb",
    params:
        dom_indir = lambda wildcards, input: dirname(input.dominant_cluster),
        script = join(ROOT_DIR, "workflow/notebooks/lineage_clones/clone_cluster_count_embed", "aggregate_clone_meta.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
        input_id = "Input" if "Input" in samples.index else "None"
    shell:
        "papermill -p dom_indir {params.dom_indir} -p clone_change_f {input.clone_change_f} -p outdir {params.outdir} -p input_id {params.input_id} {params.script} {output.note}"

# rule merge_enrich_and_lineage_and_input:
#     input:
#         clone_change_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/clone_change.csv",
#         dominant_cluster = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/dominant_clone_clust/combinedConditions_dominant_cluster.csv",
#         se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/SE.rds",
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/merge_enrich_and_lineage_and_input/out.ipynb",
#     params:
#         dom_indir = lambda wildcards, input: dirname(input.dominant_cluster),
#         script = join(ROOT_DIR, "workflow/notebooks/lineage_clones/clone_cluster_count_embed", "aggregate_clone_meta.ipynb"),
#         outdir = lambda wildcards, output: dirname(output.note)
#     shell:
#         "papermill -p se_f {input.se_f} -p dom_indir {params.dom_indir} -p enrichment_f {input.clone_change_f} -p outdir {params.outdir} {params.script} {output.note}"


# Overlay mt variants on embeddings
rule mtVarsPlot:
    input:
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/SE.rds",
        mt_cells = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/mt_clones/thr__{t}_rt__{rt}/cells_meta.tsv"
    output:
        note =  "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/thr__{t}_rt__{rt}/annotation_clones/de_clone_btwnvars_RNA_af/mtPlots.ipynb",
    params:
        mt_cells = lambda wildcards, input: input.mt_cells.replace("cells_meta.tsv", "cells_mt.tsv"),
        rscript= join(ROOT_DIR, "workflow/notebooks/lineage_clones/mtVarsPlot.ipynb"), # The script defaults to the granja data
        outdir = lambda wildcards, output: dirname(output.note),
    # test_use="wilcox",
    # latent_vars="NULL",
    threads: 8
    shell: "papermill -p se_f {input.se_f} -p cells_meta_f {input.mt_cells} -p outdir {params.outdir} {params.rscript} {output.note}"

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



# use rule dominant_clone_clust_in_clone from annCloMod as annClo_dominant_clone_clust_in_clone with:
#     params:
#         cluster_names_f = cluster_names_f
    # input:
    #     se_meta="{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/se_cells_meta.tsv",
    #     cluster_names_f = lambda wildcards: get_clust_markers
    # output:
    #     note= "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/dominant_clone_clust/dominant.ipynb",
    #     results = report(multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/dominant_clone_clust/",
    #         "dominant_cluster.csv", "cluster_clone_meta.csv", "cluster_clone_counts_normalized.csv",
    #         "combinedConditions_dominant_cluster.csv", "combinedConditions_cluster_clone_meta.csv", "combinedConditions_cluster_clone_counts_normalized.csv"))

## Use annotation_clones to run DE on clones
# use rule summary_TF_largeClones from annCloMod as annClo_summary_TF_largeClones with:
#     input:
#         se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/se.clonesfilt.rds",
#         de = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb"
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/output_Summary.ipynb"
#

# ## Use annotation_clones to run DE on clones
# use rule summary_TF_largeClones from annCloMod as knn_annClo_summary_TF_largeClones with:
#     input:
#         se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/se.clonesfilt.rds",
#         de = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb"
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/output_Summary.ipynb",


## Use annotation_clones to run DE on clones
# use rule summary_TF_largeClones from annCloMod as annClo_summary_TF_largeClones with:
#     input:
#         se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/se.clonesfilt.rds",
#         de = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb"
#     output:
#         out = multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/",
#                        "output_Summary.ipynb", "dotplot.allDonors.top.pdf", "dotplot.allDonors.clones.top.pdf",
#                        "heatmap.allDonors.split.top.pdf", "large.clones.cdf.png"),
#
# ## Use annotation_clones to run DE on clones
# use rule summary_TF_largeClones from annCloMod as knn_annClo_summary_TF_largeClones with:
#     input:
#         se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/se.clonesfilt.rds",
#         de = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb"
#     output:
#         out = report(multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/",
#                          "output_Summary.ipynb", "dotplot.allDonors.top.pdf", "dotplot.allDonors.clones.top.pdf",
#                          "heatmap.allDonors.split.top.pdf", "large.clones.cdf.png"),
#                          category="lineage and nuclear", subcategory="Method: knn Variants: {variants} resolution: {kparam} minPct: {minPct} logThresh: {logThresh}"),
#         donorDot = report(expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{{variants}}/knn/kparam_{{kparam}}/gff_{{gff}}/annotation_clones/DE_large/minPct_{{minPct}}__logThresh_{{logThresh}}/cdf_thresh__{{cdf_thresh}}/donor{d}_TF/dot.top.png",
#                                  d=range(config["N_DONORS"])),
#                           category="lineage and nuclear", subcategory="Method: KNN Variants: {variants} k resolution: {kparam} minPct: {minPct} logThresh: {logThresh}")


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

# use rule btwnClone_runGSEA_sepDonors from cloneChangeMod as cloneChange_btwnClone_runGSEA_sepDonors with:
#     input:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/de.ipynb"
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/sepDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb"
#
#
# use rule btwnClone_summaryGSEA_sepDonors from cloneChangeMod as cloneChange_btwnClone_summaryGSEA_sepDonors with:
#     input:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/sepDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb"
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/sepDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
#         #note2= "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/sepDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb"
#


############################################
## Somatic variants and clones:
############################################
params["outdir"] = config["outdir"]
params["prefix"] = config["prefix"]
#params[]
# module somvarsMod:
#     snakefile: "./rules/somatic_variants/vcf_bam_to_cells.smk"
#     config: params
#
# use rule * from somvarsMod
#/data/Mito_Trace/output/somatic_variants/CHIP_dec172021/aggregate/needle_post/peaks_chip_genes/CHIP_b1_Control/scPileupVars

# def needle_af(wildcards):
#     w = wildcards
#     samplename = f"{config['experiment']}_{w.sample}"
#     curr_dir = join(config["outdir"], "somatic_variants", config["somvars_dir"], f"aggregate/needle_post/peaks_{w.regions}/{samplename}/cells_vars/vcfpad_1")
#     return join(curr_dir, "af.pileup.tsv")
#
# def needle_ref(wildcards):
#     w = wildcards
#     samplename = f"{config['experiment']}_{w.sample}"
#     curr_dir = join(config["outdir"], "somatic_variants", config["somvars_dir"], f"aggregate/needle_post/peaks_{w.regions}/{samplename}/cells_vars/vcfpad_1")
#     return  join(curr_dir, "af.ref.pileup.tsv")
#
# def vars_anno(wildcards):
#     w = wildcards
#     return join(config["outdir"], "somatic_variants", config["somvars_dir"], f"aggregate/needle_post/peaks_{w.regions}/variants.annotate.gene.vcf")
#
#
# # use rule pileup_to_cell_vars from somvarsMod with:
# #     input:
# #     output:
# #         ref = "{config['outdir']}/aggregate/needle_post/peaks_{peaks}/{sample}/cells_vars/vcfpad_1/af.ref.pileup.tsv",
# #         af = "{out_dir}/aggregate/needle_post/peaks_{peaks}/{sample}/cells_vars/vcfpad_1/af.pileup.tsv"
#
# rule merge_som_clones:
#     input:
#         af = needle_af,
#         ref = needle_ref,
#         cells_meta = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
#         #clones = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/_clone_complete.txt",
#         #ref= "{out_dir}/aggregate/needle_post/peaks_{peaks}/{sample}/cells_vars/vcfpad_1/af.ref.pileup.tsv",
#         #af= "{out_dir}/aggregate/needle_post/peaks_{peaks}/{sample}/cells_vars/vcfpad_1/af.pileup.tsv"
#     output:
#           note="{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/needle/peaks_{regions}/{sample}/som_clones.ipynb"
#           # figs=multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/needle/peaks_{regions}/{sample}/",
#           #          "clone_altAndRef_numCells.png","clone_altAndRef_sumCounts.png", "donor_altAndRef_numCells.png","donor_altAndRef_sumCounts.png")
#     params:
#         indir = lambda wildcards, input: dirname(input.af),
#         outdir = lambda wildcards, output: dirname(output.note),
#         #cells_meta = lambda wildcards, input: join(dirname(dirname(input.clones)), "cells_meta.tsv"),
#         script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/vars_to_clones.ipynb")
#     shell: "papermill -p cells_meta_f {input.cells_meta} -p indir {params.indir} -p condition {wildcards.sample} -p outdir {params.outdir} {params.script} {output.note}"
#
#
# rule som_clones_samples:
#     input:
#          expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{{variants}}/knn/kparam_{{kparam}}/somatic_variants/needle/peaks_{{regions}}/{sample}/som_clones.ipynb",
#              sample = samples.index)
#          #note="{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/needle/peaks_{regions}/{sample}/som_clones.ipynb",
#     output:
#           vars="{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/needle/peaks_allSamples_{regions}/allSamples.variants.tsv",
#           stats="{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/needle/peaks_allSamples_{regions}/allSamples.variants.stats.tsv"
#     run:
#         allVars = []
#         for i in input:
#             curr_dir = dirname(i)
#             curr_f = join(curr_dir, "clones_alt_numCells.csv")
#             if os.path.exists(curr_f):
#                 curr_vars = pd.read_csv(curr_f)
#                 allVars.append(curr_vars)
#         if len(allVars) != 0:
#             pd.concat(allVars).to_csv(output.vars)
#         else:
#             os.system(f"touch {output.vars}")
#
# rule som_dendro:
#     input:
#         note = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{{variants}}/knn/kparam_{{kparam}}/somatic_variants/needle/peaks_{{regions}}/{sample}/som_clones.ipynb",
#              sample = samples.index),
#         vars_f = vars_anno, #"{outdir}/regions_{peaks}/gatk_mutect/post/variants.annotate.gene.vcf"
#         dendro_f = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{{variants}}/knn/kparam_{{kparam}}/barcodes/btwnClones_dendro_dt_{{dendro_thresh}}/donor{d}.mean.csv",
#                           d= np.arange(config["N_DONORS"])),
#          #"{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.mean.csv",
#         cells_meta = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
#     output:
#         note="{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/somatic_variants/needle/peaks_{regions}/som_dendro_{dendro_thresh}/som_clones.ipynb",
#     params:
#         indir = lambda wildcards, input: dirname(dirname(input.note[0])),
#         outdir = lambda wildcards, output: dirname(output.note),
#         dendro_indir = lambda wildcards, input: dirname(input.dendro_f[0]),
#         #cells_meta = lambda wildcards, input: join(dirname(dirname(input.clones)), "cells_meta.tsv"),
#         script = join(ROOT_DIR, "workflow/notebooks/somatic_variants_clones/dendro_vars_to_clones.ipynb")
#     shell: "papermill -p cells_meta_f {input.cells_meta} -p indir {params.indir} -p dendro_indir {params.dendro_indir} -p outdir {params.outdir} -p vars_f {input.vars_f} {params.script} {output.note}"


## Somatic variants merge
# Run for gatk_mutect and needle.
module somvarsMod:
    snakefile: "./rules/somatic_variants/affilt_merge_somvars_clones.smk"
    config: config


#use rule * from somvarsMod
module mtCloneMod:
    snakefile: "./rules/mt_as_clone.smk"
    config: params

use rule * from mtCloneMod as mtclone_*


module cloneShiftMod:
    snakefile: "./rules/clonal_shift.smk"
    config: params

use rule * from cloneShiftMod as cloneshift_*

module indClonesMod:
    snakefile: "./rules/individual_clones.smk"
    config: params

use rule * from indClonesMod as indclones_*

# def get_script(wildcards):
#     clone_type = wildcards.clone_type
#     if clone_type == "clones":
#         return join(ROOT_DIR, "workflow/notebooks/clonal_shifts/hypergeometric_clones.ipynb")
#     elif clone_type == "mt_den":
#         return join(ROOT_DIR, "workflow/notebooks/clonal_shifts/hypergeometric_mt_as_clones_dendro.ipynb")
#     elif clone_type == "mt_bin":
#         return join(ROOT_DIR, "workflow/notebooks/clonal_shifts/hypergeometric_mt_as_clones.ipynb")
#     elif clone_type == "dendro_bc":
#         return join(ROOT_DIR, "workflow/notebooks/clonal_shifts/dendro.ipynb")
#     raise ValueError("clone_type variable")
#     return
#
# #"{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/mt_as_clones/variants_{variants}/
# #"clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/",
# rule clonalshift_mt_as_clones:
#     input:
#         clones = "{outdir}/mt_as_clones/variants_{variants}/dendro/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/best_params_save.ipynb",
#         se_meta = "{outdir}/annotation_clones/se_cells_meta_labels.tsv"
#     output:
#         "{outdir}/clonal_shifts/mt_as_clones/results/{condition}/variants_{variants}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
#     params:
#         script = get_script,
#         outdir = lambda wildcards, output: dirname(output[0])
#     shell: "papermill -p se_cells_meta_f {input.se_meta} -p outdir {params.outdir} {params.script} {output}"
#
#
# rule clonalshift_mt_as_clones_dendro:
#     input:
#         clones = "{outdir}/mt_as_clones/variants_{variants}/dendro/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/dendro_mt_clones.ipynb",
#         se_meta = "{outdir}/annotation_clones/se_cells_meta_labels.tsv"
#     output:
#         "{outdir}/clonal_shifts/mt_as_clones_dendro/results/{condition}/variants_{variants}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
#     params:
#         indir = lambda wildcards, input: dirname(input.clones),
#         script = get_script,
#         outdir = lambda wildcards, output: dirname(output[0])
#     shell: "papermill -p indir {params.indir} -p se_cells_meta_f {input.se_meta} -p outdir {params.outdir} {params.script} {output}"
#
#
# rule clonalshift_clones:
#     input:
#         se_meta = "{outdir}/annotation_clones/variants_{variants}/se_cells_meta_labels.tsv"
#     output:
#         "{outdir}/clonal_shifts/clones/results/{condition}/variants_{variants}/output.ipynb"
#     params:
#         script = get_script,
#         outdir = lambda wildcards, output: dirname(output[0])
#     shell: "papermill -p se_cells_meta_f {input.se_meta} -p outdir {params.outdir} {params.script} {output}"
#
#
# rule clonalshift_dendro_bc:
#     input:
#         se_meta = "{outdir}/clones/annotation_clones/variants_{variants}/se_cells_meta_labels.tsv",
#         barcodes_dir = "{outdir}/clones/variants_{variants}/knn/kparam_{kparams}/barcodes/btwnClones_dendro_dt_{dt}/",
#     output:
#         "{outdir}/clones/variants_{variants}/knn/kparam_{kparams}/barcodes/btwnClones_dendro_dt_{dt}/results/{condition}/output.ipynb"
#     params:
#         script = get_script,
#         outdir = lambda wildcards, output: dirname(output[0])
#     shell: "papermill -p se_cells_meta_f {input.se_meta} -p outdir {params.outdir} {params.script} {output}"


# use rule nuclear_and_mtclone_counts from mtcloneMod as mtclone_nuclear_and_mtclone_counts with:
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_{variant}/anno_mt_clones_thresh/anno_mt_clones.ipynb",
#
#
# use rule nuclear_and_mtclone_counts from mtcloneMod as mtclone_ with:
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_A2_black/mt_as_clones/out.ipynb"

#
# use rule finalize from mtCloneMod as mtclone_finalize with:
#     output:
#         out = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/mt_as_clones/variants_{variants}/.complete.txt"
