from os.path import join, dirname
from src.utils.data_io import sparse_to_cellranger_fragments
import pandas as pd
from src.config import ROOT_DIR
import os
from src.utils.parse_config import read_config_file
#configfile: "/data2/mito_lineage/parameters/DUPI_april08_2021/mttrace_mtasnucl.yaml"
#params = read_config_file(config["config"])


# Merge config + params together, with config the default params used
# for p in params:
#     if p not in config:
#         config[p] = params[p]

# cfg_anno = config['annotations']
# extrnl = cfg_anno["name"]
#samples = pd.read_table(join(ROOT_DIR, config["samples_meta"]), dtype=str,sep=',').set_index(["sample_name"], drop=False)
#res = config["results"]

#
# rule all:
#     input:
#         expand("{outdir}/annotation_clones/SE.rds",
#                outdir=config['outdir'], prefix=config["prefix"]),
#         expand("{outdir}/annotation_clones/DE/donor{n}/DE.ipynb",
#                 outdir=config['outdir'], prefix=config["prefix"],
#                 n=config['N_DONORS']),
#         expand("{outdir}/annotation_clones/DE_TF/donor{n}/DE_TF.ipynb",
#                 outdir=config['outdir'], prefix=config["prefix"],
#                 n=config['N_DONORS']),
#
#         expand("{outdir}/annotation_clones/DE/donor{n}/GSEA/clusters/GSEA.ipynb",
#                 outdir=config['outdir'], prefix=config["prefix"],
#                 n=config['N_DONORS']),
#
#         expand("{outdir}/annotation_clones/allDonors/DE/DE.ipynb",
#                 outdir=config['outdir'], prefix=config["prefix"],
#                 ),


def get_anno_integrate(wildcards):
    print(join(config["anno_res"], "mergedSamples","allSamples.integrated.rds"))
    return join(config["anno_res"], "mergedSamples","allSamples.integrated.rds")


rule addClones:
    input:
        noc = get_anno_integrate, #"{outdir}/annotation_clones/mergedSamples/allSamples.integrated.rds", "{anno_indir}/mergedSamples/allSamples.integrated.rds", #
        clones = "{outdir}/cells_meta.tsv",#"{outdir}/pipeline/published/{prefix}/data/clones/clones.txt",
    output:
        se_f = "{outdir}/annotation_clones/SE.rds",
        note = "{outdir}/annotation_clones/addClones.ipynb"
    params:
        outdir = lambda wildcards, output: dirname(output[0]),
        rscript= join(ROOT_DIR, "R_scripts/annotation_clones/addClones.v01.vCurrent.ipynb"), # The script defaults to the granja data
    shell: "papermill -p cells_meta_f {input.clones} -p se_f {input.noc} -p outdir {params.outdir} {params.rscript} {output.note}"


rule plotMarkers:
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
    output:
        note = "{outdir}/annotation_clones/markers/markers.ipynb",
    params:
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/umap_immune_markers.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
        markers_f = config["markers_f"]
    shell: "papermill -p se_f {input.se_f} -p outdir {params.outdir} -p markers_f {params.markers_f} {params.rscript} {output.note}"

rule counts_clones:
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
    output:
        multiext("{outdir}/annotation_clones/clone_counts/minCellConds_{min_cells}/",
            "clone_sizes.ipynb", "clone_sizes.csv", "clone_sizes_norm.csv")
    params:
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/counts_inClones.ipynb"),
        outdir = lambda wildcards, output: dirname(output[0]),
        sample_names = ",".join(config["samples"].index),
    shell:
        "papermill -p se_f {input.se_f} -p outdir {params.outdir} -p sample_names {params.sample_names} -p minCell {wildcards.min_cells} {params.rscript} {output[0]}"


rule runDE_enrich:
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
        enrich_f = expand("{{outdir}}/enrichmentNorm_donor{d}.csv",
                          d=range(config["N_DONORS"])),
    output:
        note="{outdir}/annotation_clones/DE/DE.ipynb"
    params:
        #d = lambda wildcards: wildcards.d,
        enrich_d = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        N_DONORS = config["N_DONORS"],
        rscript= join(ROOT_DIR, "R_scripts/annotation_clones/addEnrichment.vCurrent.ipynb")
    shell: "papermill -p enrich_d {params.enrich_d} -p se_f {input.se_f} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} {params.rscript} {output.note}"


########################################################################################################################
## TF analysis on larger clones
########################################################################################################################
rule setup_motifs_largeClones:
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
    output:
        note="{outdir}/annotation_clones/DE_large/output.ipynb",
        largeclones_se = "{outdir}/annotation_clones/DE_large/se.clonesfilt.rds"
    threads: 16
    params:
        outdir = lambda wildcards, output: dirname(output[0]),
        n_donors = config["N_DONORS"],
        rscript= join(ROOT_DIR, "R_scripts/annotation_clones/setup_motifs_large_clones.ipynb")
    shell: "papermill -p se_f {input.se_f} -p outdir {params.outdir} -p n_donors {params.n_donors} {params.rscript} {output.note}"


rule runDE_TF_largeClones:
    input:
        se_f = "{outdir}/annotation_clones/DE_large/se.clonesfilt.rds"
    output:
        note = "{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb",
        out2 = report(directory(expand("{{outdir}}/annotation_clones/DE_large/minPct_{{minPct}}__logThresh_{{logThresh}}/donor{d}_TF",
                                        d = range(config["N_DONORS"]))),
                      patterns=["*png", "*csv"],
                      category="lineage", subcategory="donor {d}"),
    threads: 8
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        n_donors = config["N_DONORS"],
        rscript= join(ROOT_DIR, "R_scripts/annotation_clones/DE_clones.TF.vCurrent.ipynb"),
        #min_pct = lambda wildcards
    shell: "papermill -p indir {params.indir} -p outdir {params.outdir} -p min_pct {wildcards.minPct} -p logfc_thresh {wildcards.logThresh} -p n_donors {params.n_donors} {params.rscript} {output.note}"


rule summary_TF_largeClones:
    input:
        de = "{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb",
        se_f = "{outdir}/annotation_clones/DE_large/se.clonesfilt.rds"
    output:
        out = multiext("{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/",
                      "output_Summary.ipynb", "dotplot.allDonors.top.pdf", "dotplot.allDonors.clones.top.pdf",
                      "heatmap.allDonors.top.pdf",
                      "heatmap.allDonors.split.top.pdf", "large.clones.cdf.png"),
        donorDot = expand("{{outdir}}/annotation_clones/DE_large/minPct_{{minPct}}__logThresh_{{logThresh}}/cdf_thresh__{{cdf_thresh}}/donor{d}_TF/dot.top.png",
                        d=range(config["N_DONORS"])),
        #note="{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/output_Summary.ipynb",
    params:
        indir = lambda wildcards, input: dirname(input.de),
        se_indir = lambda wildcards, input: dirname(input.se_f),
        n_donors = config["N_DONORS"],
        cdf_thresh = lambda wildcards: wildcards.cdf_thresh,
        rscript= join(ROOT_DIR, "R_scripts/annotation_clones/DE_clones.TF.summaryPlot.ipynb")
    shell: "papermill -p indir {params.indir} -p se_indir {params.se_indir} -p n_donors {params.n_donors} -p cdf_thresh {params.cdf_thresh} {params.rscript} {output[0]}" #{output[0]}

rule output_largeclones:
    """ Summarizes TF results across all clonal comparisons. 
    """
    input:
        de = "{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb",
        se_f = "{outdir}/annotation_clones/DE_large/se.clonesfilt.rds"
    output:
        note = "{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/TFactivity.ipynb",
        out = multiext("{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/",
                       "dotplot.allDonors.clones.topPvals.pdf", "dotplot.allDonors.clones.topPvals.noZ.pdf",
                       "clone.TFactivity.topPvals.txt")
    params:
        indir = lambda wildcards, input: dirname(input.de),
        se_indir = lambda wildcards, input: dirname(input.se_f),
        n_donors = config["N_DONORS"],
        cdf_thresh = lambda wildcards: wildcards.cdf_thresh,
        rscript = join(ROOT_DIR,"R_scripts/annotation_clones/DE_clones.TF.summaryPlot.DonorConcat.Development.ipynb")
    shell: "papermill -p indir {params.indir} -p se_indir {params.se_indir} -p n_donors {params.n_donors} -p cdf_thresh {params.cdf_thresh} {params.rscript} {output.note}"

rule runGSEA:
    """ GSEA on large clones DE
            
        TODO: implement
    """
    input:
        DE_out_path = "{outdir}/annotation_clones/DE/donor{n}/DE.ipynb"
    output:
        "{outdir}/annotation_clones/DE/donor{n}/GSEA/clusters/GSEA.ipynb",
    params:
        input = lambda wildcards, input: join(dirname(input[0]), "clusters"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "R_scripts/annotations/clones/runGSEA_clones.ipynb")
    shell: "papermill -p DE.out.path {params.input} -p export.path {params.output} {params.rscript} {output[0]}"


# rule runDE_TF:
#     input:
#         integrated = "{outdir}/annotation_clones/{prefix}/mergedSamples/allSamples.integrated.rds",
#     output:
#         "{outdir}/annotation_clones/{prefix}/mergedSamples/DE_TF/DE_TF.ipynb",
#     params:
#         outdir = lambda wildcards, output: dirname(output[0]),
#         rscript= join(ROOT_DIR, "R_scripts/annotations/DE_TF.ipynb"), # The script defaults to the granja data
#         sample_names = ",".join(samples.index),
#         comps_f = get_comps()
#     shell: "papermill -p integrated_f {input} -p outdir {params.outdir} -p sample_names {params.sample_names} -p comps_f {params.comps_f:q} {params.rscript} {output[0]}"
#

# rule overlay_cells_meta:
#     input:
#         "{outdir}/annotation_clones/{prefix}/allSamples.integrated.rds",
#     output:
#         "{outdir}/annotation_clones/{prefix}/donors.png",
#         "{outdir}/annotation_clones/{prefix}/clones.png"
#     params:
#         indir = lambda wildcards, input: dirname(input[0]),
#         outdir = lambda wildcards, output: dirname(output[0]),
#         cells_meta = config["cells_meta"],
#         rscript= join(ROOT_DIR, "R_scripts/annotations/lineageNuclear.ipynb"), # The script defaults to the granja data
#
#     shell: "papermill -p cellr_in {params.indir} -p outdir {params.outdir} -p cells_meta {params.cells_meta} {params.rscript} {output[0]}"
#

## Till Here
########################################################################################################################

########################################################################################################################
# Run DE for genes and peak and GSEA (for genes) in different comparisons.
# Comparisons include either across nuclear clusters, across conditions within each cluster,
# and across conditions within each clone.
########################################################################################################################

## RUN DE for gene activity and peaks as markers.
## Use the demultiplexed cells from addClones.
def get_btwnClust_rscript(wildcards):
    if wildcards.assay == "RNA":
        return join(ROOT_DIR, "R_scripts/annotation_clones/DE_genes_btwnClusters.ipynb")
    else:
        return join(ROOT_DIR, "R_scripts/annotation_clones/DE_peaks_btwnClusters.ipynb")
    return

rule btwnClust_DE:
    """ Compares clusters to each other. For now works with Gene Activity
    
    Not tested with peaks as gene activity yet.
    """
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
        rscript = get_btwnClust_rscript
    output:
        note =  "{outdir}/annotation_clones/de_btwnclust_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}/pthresh_{p_thresh}.ipynb",
    params:
        outdir = lambda wildcards, output: dirname(output.note),
        assay=lambda wildcards: wildcards.assay,
        minPct=lambda wildcards: wildcards.btwnMinpct,
        logfcthresh= lambda wildcards: wildcards.logfc_threshold,
        top_de=3,
        samples = ",".join(config["samples"].index),
        p_thresh=lambda wildcards: wildcards.p_thresh,
    shell: "papermill -p se_f {input.se_f} -p outdir {params.outdir} -p top_de {params.top_de} -p sample_names {params.samples}  {input.rscript} {output.note}"

rule runGSEA_btwnClust:
    input:
        DE_out_path = "{outdir}/annotation_clones/de_btwnclust_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}/pthresh_{p_thresh}.ipynb",
    output:
        note= "{outdir}/annotation_clones/de_btwnclust_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}/de_p{p_thresh}_GSEA_pthresh_{gsea_pval}/GSEA.ipynb",
    #note="{outdir}/data/annotation/gff_{gff}/mergedSamples/DE/GSEA/clusters/GSEA.ipynb",
    params:
        input = lambda wildcards, input: join(dirname(input[0]), "btwnClust"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/runGSEA_btwnClust.ipynb"),
        gsea_dir = join(ROOT_DIR, "software/Bioinformatics_Tools/"), #/ is necessary
        gsea_pval = lambda wildcards: wildcards.gsea_pval
        #conda_env = "../envs/gsea_manual.yml"
    conda:
        "../envs/gsea_manual.yml" #environment with clusterprofiler
    shell: "papermill -p DE.out.path {params.input} -p gsea_pthresh {params.gsea_pval} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"

rule summaryGSEA_btwnClust:
    input:
        #"{outdir}/annotation_clones/de_btwnclust_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_p{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
        note = "{outdir}/annotation_clones/de_btwnclust_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}/de_p{p_thresh}_GSEA_pthresh_{gsea_pval}/GSEA.ipynb",
    output:
        note= "{outdir}/annotation_clones/de_btwnclust_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}/de_p{p_thresh}_GSEA_pthresh_{gsea_pval}/summary.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]), #join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/summarizeGSEA.ipynb"),
    shell: "papermill -p export_path {params.input} {params.rscript} {output.note}"


## RUN DE for gene activity and peaks as markers.
## Use the demultiplexed cells from addClones.
def get_btwnCond_rscript(wildcards):
    if wildcards.assay == "RNA":
        return join(ROOT_DIR, "R_scripts/annotation_clones/DE_genes_btwnConds.ipynb")
    else:
        return join(ROOT_DIR, "R_scripts/annotation_clones/DE_peaks_btwnConds.ipynb")
    return

rule btwnCond_DE:
    """ Compares clusters to each other. For now works with Gene Activity
    
    Not tested with peaks as gene activity yet.
    """
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
    output:
        note =  "{outdir}/annotation_clones/de_btwncond_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}/pthresh_{p_thresh}.ipynb",
    params:
        rscript= join(ROOT_DIR, "R_scripts/annotation_clones/DE_genes_btwnConditions.ipynb"), # The script defaults to the granja data
        samples = ",".join(config["samples"].index),
        outdir = lambda wildcards, output: f"{dirname(output.note)}_pthresh_{wildcards.p_thresh}", #"{outdir}/annotation_clones/de_btwnConds_{assay}/minPct_{btwnMinpct}_logfc{logfcthresh}/pthresh_{p_thresh}",
        assay=lambda wildcards: wildcards.assay,
        minPct=lambda wildcards: wildcards.btwnMinpct,
        logfcthresh= lambda wildcards: wildcards.logfc_threshold,
        top_de=3,
        p_thresh=lambda wildcards: wildcards.p_thresh,
        # test_use="wilcox",
        # latent_vars="NULL",
    threads: 8
    shell: "papermill -p se_f {input.se_f} -p p_thresh {params.p_thresh} -p outdir {params.outdir} -p top_de {params.top_de} -p sample_names {params.samples} {params.rscript} {output.note}"


rule runGSEA_btwnCond:
    input:
        DE_out_path = "{outdir}/annotation_clones/de_btwncond_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}/pthresh_{p_thresh}.ipynb",
    output:
        note= "{outdir}/annotation_clones/de_btwncond_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}_p{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
    params:
        input = lambda wildcards, input: join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/runGSEA_btwnCond.ipynb"),
        gsea_dir = join(ROOT_DIR, "software/Bioinformatics_Tools/"), #/ is necessary
        gsea_pval = lambda wildcards: wildcards.gsea_pval,
        stat_col = lambda wildcards: wildcards.stat_col,
        prefilter = lambda wildcards: wildcards.prefilter,
        padjmethod = lambda wildcards: wildcards.padjmethod,
        pthresh = lambda wildcards: wildcards.p_thresh
    #conda_env = "../envs/gsea_manual.yml"
    conda:
        "../envs/gsea_manual.yml" #environment with clusterprofiler
    shell: "papermill -p padjmethod {params.padjmethod} -p pthresh {params.pthresh} -p gsea_pthresh {params.gsea_pval} -p prefilter {params.prefilter} -p stat_col {params.stat_col:q} -p DE.out.path {params.input} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"

rule summaryGSEA_btwnCond:
    input:
        "{outdir}/annotation_clones/de_btwncond_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}_p{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
    output:
        note="{outdir}/annotation_clones/de_btwncond_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}_p{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]), #join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/summarizeGSEA.ipynb"),
    shell: "papermill -p export_path {params.input} {params.rscript} {output.note}"


rule inClone_btwnCond_DE:
    """ Compares clusters to each other. For now works with Gene Activity
    
    Not tested with peaks as gene activity yet.
    """
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
    output:
        note =  "{outdir}/annotation_clones/de_clone_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/de.ipynb",
    params:
        rscript= join(ROOT_DIR, "R_scripts/annotation_clones/DE_genes_inClones_btwnConditions.ipynb"), # The script defaults to the granja data
        samples = ",".join(config["samples"].index),
        outdir = lambda wildcards, output: dirname(output.note), #"{outdir}/annotation_clones/de_btwnConds_{assay}/minPct_{btwnMinpct}_logfc{logfcthresh}/pthresh_{p_thresh}",
        minPct=lambda wildcards: wildcards.btwnMinpct,
        logfcthresh= lambda wildcards: wildcards.logfc_threshold,
        top_de=3,
        p_thresh=lambda wildcards: wildcards.p_thresh,
    # test_use="wilcox",
    # latent_vars="NULL",
    threads: 8
    shell: "papermill -p se_f {input.se_f} -p min_cells {wildcards.min_cells} -p p_thresh {params.p_thresh} -p outdir {params.outdir} -p top_de {params.top_de} -p sample_names {params.samples} {params.rscript} {output.note}"

rule runGSEA_inClone_btwnCond:
    input:
        DE_out_path = "{outdir}/annotation_clones/de_clone_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/de.ipynb",
    output:
        note= "{outdir}/annotation_clones/de_clone_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]), #join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/runGSEA_inClones_btwnConditions.ipynb"),
        gsea_dir = join(ROOT_DIR, "software/Bioinformatics_Tools/"), #/ is necessary
        gsea_pval = lambda wildcards: wildcards.gsea_pval,
        stat_col = lambda wildcards: wildcards.stat_col,
        prefilter = lambda wildcards: wildcards.prefilter,
        padjmethod = lambda wildcards: wildcards.padjmethod,
        pthresh = lambda wildcards: wildcards.p_thresh,
        pattern = "*clone_.*counts.*"
    #conda_env = "../envs/gsea_manual.yml"
    conda:
        "../envs/gsea_manual.yml" #environment with clusterprofiler
    shell: "papermill  -p DE.out.path {params.input} -p export.path {params.output} -p pattern {params.pattern:q} -p padjmethod {params.padjmethod} -p pthresh {params.pthresh} -p gsea_pthresh {params.gsea_pval} -p prefilter {params.prefilter} -p stat_col {params.stat_col:q} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"

rule summaryGSEA_inClone_btwnCond:
    input:
        "{outdir}/annotation_clones/de_clone_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
    output:
        note="{outdir}/annotation_clones/de_clone_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]), #join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        script = join(ROOT_DIR, "R_scripts/annotation_clones/summarizeGSEA.ipynb"),
    shell: "papermill -p export_path {params.input} {params.script} {output.note}"


####################################
## Clones-clusters analysis by size.
####################################

rule se_meta:
    """Prepare clone-by-cluster counts for umap and hypergeometric test"""
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
    output:
        se_meta = "{outdir}/annotation_clones/se_cells_meta.tsv",
        note = "{outdir}/annotation_clones/se_cells_meta.ipynb",
    params:
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/clusters_cells_meta.ipynb"),
    shell: "papermill -p se_f {input.se_f} {params.rscript} {output.note}"


rule dominant_clone_clust_in_clone:
    input:
        se_meta = "{outdir}/annotation_clones/se_cells_meta.tsv",
    output:
        note= "{outdir}/annotation_clones/dominant_clone_clust/dominant.ipynb",
        results = report(multiext("{outdir}/annotation_clones/dominant_clone_clust/",
            "dominant_cluster.csv", "cluster_clone_meta.csv", "cluster_clone_counts_normalized.csv",
            "combinedConditions_dominant_cluster.csv", "combinedConditions_cluster_clone_meta.csv", "combinedConditions_cluster_clone_counts_normalized.csv"))
    params:
        outdir = lambda wildcards, output: dirname(output.note),
        script = join(ROOT_DIR, "src", "clone_cluster_count_embed", "dominant_clust_in_clone.ipynb"),
        samples = ",".join(config["samples"].index),
        cluster_names_f= config.get("cluster_names_f", "None")
    shell:
        "papermill -p outdir {params.outdir} -p se_cells_meta_f {input.se_meta} -p sample_names {params.samples} -p cluster_names_f {params.cluster_names_f:q} {params.script} {output.note}"

rule count_clust_embed:
    input:
        dominant = "{outdir}/annotation_clones/dominant_clone_clust/dominant.ipynb"
    output:
        note = "{outdir}/annotation_clones/clone_clust_embed/tsne_perp{perp}_donperp{donperp}/embed.ipynb",
        sepDonors = "{outdir}/annotation_clones/clone_clust_embed/tsne_perp{perp}_donperp{donperp}/combinedConditions_combinedDonors/umap.pdf",
        combinedDonors = "{outdir}/annotation_clones/clone_clust_embed/tsne_perp{perp}_donperp{donperp}/combinedDonors/umap.pdf"
        #"{outdir}/annotation_clones/clone_clust_umap/"
    params:
        indir = lambda wildcards, input: dirname(input.dominant),
        outdir = lambda wildcards, output: dirname(output.note),
        samples = ",".join(config["samples"].index),
        script = join(ROOT_DIR, "src/clone_cluster_count_embed/run_count_clust_embed.ipynb"),
    shell:  "papermill -p indir {params.indir} -p outdir {params.outdir} -p sample_names {params.samples} -p perplexity {wildcards.perp} {params.script} {output.note}"


rule cluster_clone_hypergeom:
    """Hypergeometric distribution for clones and clusters"""
    input:
        se_meta = "{outdir}/annotation_clones/se_cells_meta.tsv",
    output:
        result="{outdir}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/result.csv",
        note="{outdir}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/result.ipynb"
    params:
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/clones_clusters_hypergeometric.ipynb")
    shell: "papermill -p out_f {output.result} -p se_cells_meta_f {input.se_meta} -p min_clone_size {wildcards.hyperMinCl} -p conds_sep {wildcards.bothConds} -p p_thresh {wildcards.pthresh} {params.rscript} {output.note}"


# # DE btwnVars
# rule mtVarsPlot:
#     input:
#         se_f = "{outdir}/annotation_clones/SE.rds",
#         mt_cells = "cells_meta.tsv" if "cells_meta" not in config else config["cells_meta"]
#     output:
#         note =  "{outdir}/annotation_clones/de_clone_btwnvars_RNA_af/mtPlots.ipynb",
#     params:
#         mt_cells = lambda wildcards, input: input.mt_cells.replace("cells_meta.tsv", "cells_mt.tsv"),
#         rscript= join(ROOT_DIR, "R_scripts/annotation_clones/mtVarsPlot.ipynb"), # The script defaults to the granja data
#         outdir = lambda wildcards, output: dirname(output.note),
#     # test_use="wilcox",
#     # latent_vars="NULL",
#     threads: 8
#     shell: "papermill -p se_f {input.se_f} -p cells_meta_f {input.mt_cells} -p outdir {params.outdir} {params.rscript} {output.note}"

rule inVar_btwnCond_DE:
    """ Compares clusters to each other. For now works with Gene Activity
    
    Not tested with peaks as gene activity yet.
    """
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
        mt_cells = "cells_meta.tsv" if "cells_meta" not in config else config["cells_meta"]
    output:
        note =  "{outdir}/annotation_clones/de_clone_btwnvars_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/de.ipynb",
    params:
        rscript= join(ROOT_DIR, "R_scripts/annotation_clones/DE_genes_inVariants_btwnConditions.ipynb"), # The script defaults to the granja data
        samples = ",".join(config["samples"].index),
        outdir = lambda wildcards, output: dirname(output.note), #f"{dirname(output.note)}_pthresh_{wildcards.p_thresh}", #"{outdir}/annotation_clones/de_btwnConds_{assay}/minPct_{btwnMinpct}_logfc{logfcthresh}/pthresh_{p_thresh}",
        minPct=lambda wildcards: wildcards.btwnMinpct,
        logfcthresh= lambda wildcards: wildcards.logfc_threshold,
        top_de=3,
        p_thresh=lambda wildcards: wildcards.p_thresh,
    # test_use="wilcox",
    # latent_vars="NULL",
    threads: 8
    shell: "papermill -p se_f {input.se_f} -p cells_meta_f {input.mt_cells} -p outdir {params.outdir} -p min_cells {wildcards.min_cells} -p p_thresh {params.p_thresh} -p top_de {params.top_de} -p sample_names {params.samples} {params.rscript} {output.note}"


rule runGSEA_inVars_btwnCond:
    input:
        DE_out_path = "{outdir}/annotation_clones/de_clone_btwnvars_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/de.ipynb",
    output:
        note= "{outdir}/annotation_clones/de_clone_btwnvars_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]), #join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/runGSEA_inClones_btwnConditions.ipynb"),
        gsea_dir = join(ROOT_DIR, "software/Bioinformatics_Tools/"), #/ is necessary
        gsea_pval = lambda wildcards: wildcards.gsea_pval,
        stat_col = lambda wildcards: wildcards.stat_col,
        prefilter = lambda wildcards: wildcards.prefilter,
        padjmethod = lambda wildcards: wildcards.padjmethod,
        pthresh = lambda wildcards: wildcards.p_thresh,
        pattern = "*clone_*"
    #conda_env = "../envs/gsea_manual.yml"
    conda:
        "../envs/gsea_manual.yml" #environment with clusterprofiler
    shell: "papermill -p pattern {params.pattern:q} -p padjmethod {params.padjmethod} -p pthresh {params.pthresh} -p gsea_pthresh {params.gsea_pval} -p prefilter {params.prefilter} -p stat_col {params.stat_col:q} -p DE.out.path {params.input} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"

rule summaryGSEA_inVars_btwnCond:
    input:
        "{outdir}/annotation_clones/de_clone_btwnvars_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
    output:
        note= "{outdir}/annotation_clones/de_clone_btwnvars_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]), #join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/summarizeGSEA.ipynb"),
    shell: "papermill -p export_path {params.input} {params.rscript} {output.note}"

