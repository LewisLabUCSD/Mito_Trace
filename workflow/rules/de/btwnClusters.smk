from os.path import join, dirname
from src.config import ROOT_DIR

########################################################################################################################
# Run DE for genes and peak and GSEA (for genes) in between clusters.
# Comparisons include either across nuclear clusters, across conditions within each cluster,
# and across conditions within each clone.
########################################################################################################################

## RUN DE for gene activity and peaks as markers.
## Use the demultiplexed cells from addClones.
def get_btwnClust_rscript(wildcards):
    if wildcards.assay == "RNA":
        return join(ROOT_DIR, "workflow/notebooks/nuclear_de/de_clusters/DE_genes_btwnClusters.ipynb")
    else:
        return join(ROOT_DIR, "workflow/notebooks/nuclear_de/de_clusters/DE_peaks_btwnClusters.ipynb")
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
        #minPct=lambda wildcards: wildcards.btwnMinpct,
        #logfcthresh= lambda wildcards: wildcards.logfc_threshold,
        top_de=3,
        samples = ",".join(config["samples"].index),
        to_donors = "TRUE"
        #p_thresh=lambda wildcards: wildcards.p_thresh,

    shell: "papermill -p se_f {input.se_f} -p outdir {params.outdir} -p top_de {params.top_de} -p minPct {wildcards.btwnMinpct} -p logfcthresh {wildcards.logfc_threshold} -p p_thresh {wildcards.p_thresh} -p sample_names {params.samples} -p to_donors {params.to_donors} {input.rscript} {output.note}"


rule runGSEA_btwnClust:
    input:
        DE_out_path = "{outdir}/annotation_clones/de_btwnclust_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}/pthresh_{p_thresh}.ipynb",
    output:
        note= "{outdir}/annotation_clones/de_btwnclust_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}/de_p{p_thresh}_GSEA_pthresh_{gsea_pval}/GSEA.ipynb",
    #note="{outdir}/data/annotation/gff_{gff}/mergedSamples/DE/GSEA/clusters/GSEA.ipynb",
    params:
        input = lambda wildcards, input: join(dirname(input[0]), "btwnClust"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "workflow/notebooks/nuclear_de/de_clusters/runGSEA_btwnClust.ipynb"),
        gsea_dir = join(ROOT_DIR, "software/Bioinformatics_Tools/"), #/ is necessary
        gsea_pval = lambda wildcards: wildcards.gsea_pval
    #conda_env = "../envs/gsea_manual.yml"
    conda:
        "../../envs/gsea_manual.yml" #environment with clusterprofiler
    shell: "papermill -p DE.out.path {params.input} -p gsea_pthresh {params.gsea_pval} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"


rule summaryGSEA_btwnClust:
    input:
        #"{outdir}/annotation_clones/de_btwnclust_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_p{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
        note = "{outdir}/annotation_clones/de_btwnclust_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}/de_p{p_thresh}_GSEA_pthresh_{gsea_pval}/GSEA.ipynb",
    output:
        note= "{outdir}/annotation_clones/de_btwnclust_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}/de_p{p_thresh}_GSEA_pthresh_{gsea_pval}/summary.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]), #join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        rscript = join(ROOT_DIR, "workflow/notebooks/utils/summarizeGSEA.ipynb"),
    shell: "papermill -p export_path {params.input} {params.rscript} {output.note}"
