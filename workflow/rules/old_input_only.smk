from os.path import join, dirname
from src.utils.data_io import sparse_to_cellranger_fragments
import pandas as pd
from src.config import ROOT_DIR
import os
from src.utils.parse_config import read_config_file



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



#####
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
