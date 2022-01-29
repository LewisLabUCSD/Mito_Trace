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
        markers_f="/data/Mito_Trace/data/processed/pbmc_markers.txt"
    shell: "papermill -p se_f {input.se_f} -p outdir {params.outdir} -p markers_f {params.markers_f} {params.rscript} {output.note}"


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
    input:
        DE_out_path = "{outdir}/annotation_clones/DE/donor{n}/DE.ipynb"
    output:
        "{outdir}/annotation_clones/DE/donor{n}/GSEA/clusters/GSEA.ipynb",
    params:
        input = lambda wildcards, input: join(dirname(input[0]), "clusters"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "R_scripts/annotations/clones/runGSEA_clones.ipynb")
    shell: "papermill -p DE.out.path {params.input} -p export.path {params.output} {params.rscript} {output[0]}"


#rule runDE_TF_Donors:


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


################
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
        #rscript = join(ROOT_DIR, "R_scripts/annotation_clones/DE_genes_btwnClusters.ipynb"), # The script defaults to the granja data
        outdir = lambda wildcards, output: dirname(output.note),
        assay=lambda wildcards: wildcards.assay,
        minPct=lambda wildcards: wildcards.btwnMinpct,
        logfcthresh= lambda wildcards: wildcards.logfc_threshold,
        top_de=3,
        samples = ",".join(config["samples"].index),
        p_thresh=lambda wildcards: wildcards.p_thresh,
        # test_use="wilcox",
        # latent_vars="NULL",
    shell: "papermill -p se_f {input.se_f} -p outdir {params.outdir} -p top_de {params.top_de} -p sample_names {params.samples}  {input.rscript} {output.note}"

rule runGSEA_btwnClust:
    input:
        DE_out_path = "{outdir}/annotation_clones/de_btwnclust_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}/pthresh_{p_thresh}.ipynb",
    #DE_out_path = "{outdir}/data/annotation/gff_{gff}/mergedSamples/DE/DE.ipynb",
    output:
        note= "{outdir}/annotation_clones/de_btwnclust_{assay}/minPct_{btwnMinpct}_logfc{logfc_threshold}_p{p_thresh}/GSEA_pthresh_{gsea_pval}/GSEA.ipynb",
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
    # shell:
    #     """
    #     conda init bash
    #     conda activate {params.conda_env}
    #     papermill -p DE.out.path {params.input} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}
    #     """
    shell: "papermill -p DE.out.path {params.input} -p pthresh {params.gsea_pval} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"


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
    # shell:
    #     """
    #     conda init bash
    #     conda activate {params.conda_env}
    #     papermill -p DE.out.path {params.input} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}
    #     """
    shell: "papermill -p padjmethod {params.padjmethod} -p pthresh {params.pthresh} -p gsea_pthresh {params.gsea_pval} -p prefilter {params.prefilter} -p stat_col {params.stat_col:q} -p DE.out.path {params.input} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"

