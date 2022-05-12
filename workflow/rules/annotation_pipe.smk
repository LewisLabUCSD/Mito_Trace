#configfile: "parameters/Zhu_Single_Cell_10X_genomics_Human2_002.yaml"
from os.path import join, dirname
from src.utils.data_io import sparse_to_cellranger_fragments
import pandas as pd
from src.config import ROOT_DIR
import os
from src.utils.parse_config import read_config_file
#configfile: "/data2/mito_lineage/parameters/DUPI_april08_2021/mttrace_mtasnucl.yaml"
params = read_config_file(config["config"])


# Merge config + params together, with config the default params used
for p in params:
    if p not in config:
        config[p] = params[p]

cfg_anno = config['annotations']
extrnl = cfg_anno["name"]
samples = pd.read_table(join(ROOT_DIR, config["samples_meta"]), dtype=str,sep=',').set_index(["sample_name"], drop=False)
#res = config["results"]
gff = params["genome_path"][config["genome"]]["gff"]


def get_cellr_input():
    out_d = cfg_anno['datadir']
    mtx_f = join(out_d, f"{cfg_anno['ID']}.mtx")
    peaks_f = join(out_d, f"{cfg_anno['ID']}.peaks.bed")
    cells_f = join(out_d, f"{cfg_anno['ID']}.cell_barcodes.txt")
    return out_d, mtx_f, peaks_f,  cells_f


rule peaks_to_fragments:
    input: get_cellr_input()
    output: "{outdir}/annotation/data/{extrnl}/{extrnl}.fragments.sort.tsv.gz"
    run: sparse_to_cellranger_fragments(input[0], input[1], input[2], out_f=output[0], n_cpus=24)


# rule sort_peaks:
#     input:
#         "{outdir}/annotation/data/{extrnl}/{extrnl}.fragments.tsv"
#     output:
#         "{outdir}/annotation/data/{extrnl}/{extrnl}.fragments.sort.tsv"
#     #shell: ""

def get_cellr_dir(wildcards):
    return join(config["mtscATAC_OUTDIR"], samples.loc[wildcards.sample, "cellr_ID"], "outs")


def get_extrnl_frags(wildcards):
    return f"{wildcards.outdir}/annotation/data/{cfg_anno['name']}/{cfg_anno['name']}.fragments.sort.tsv.gz"



#####################
## V02
#####################
rule createExpSignac:
    """ Creates a merged R SummarizedExperiment object with the external data and reference combined.
    The peaks are reduced to the overlapping sets, and if a fragments file is provided,
    Will re-count peaks after creating overlapping peak set (and changing the lengths if the overlaps extend)
    """
    output:
        "{outdir}/data/annotation/gff_{gff}/mergedSamples/createExpSignac.ipynb",
        "{outdir}/data/annotation/gff_{gff}/mergedSamples/allSamples.rds",
        report("{outdir}/data/annotation/gff_{gff}/mergedSamples/QC_01.png", category="Nuclear ATAC",
                subcategory="QC"),
        "{outdir}/data/annotation/gff_{gff}/mergedSamples/QC_01.pdf"
    params:
        indir = config["mtscATAC_OUTDIR"],
        outdir =lambda wildcards, output: dirname(output[0]),
        rscript= join(ROOT_DIR, "R_scripts/annotations/samplesCreateSignac.ipynb"),
        sample_names = ",".join(samples.index),
        samples = ",".join(samples["cellr_ID"].values),
        gff = gff
        #workdir = os.getcwd(),
    shell: "papermill -p cellr_in {params.indir} -p outdir {params.outdir} -p samples {params.samples} -p sample_names {params.sample_names} -p gff_id {params.gff} {params.rscript} {output[0]}"


def get_integrate_or_single_rscript(wildcards):
    if len(samples.index)==1:
        return join(ROOT_DIR, "R_scripts/annotations/integrateSignac_singleSample.ipynb"), # The script defaults to the granja data
    return join(ROOT_DIR, "R_scripts/annotations/integrateSignac.ipynb"), # The script defaults to the granja data



rule integrateSignac:
    input:
        a="{outdir}/data/annotation/gff_{gff}/mergedSamples/allSamples.rds",
    output:
        a="{outdir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.ipynb",
        b="{outdir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds",
        c=report(expand("{{outdir}}/data/annotation/gff_{{gff}}/mergedSamples/{f}",
                      f=["QC_02.png",
                         "integrated.merged.compare.png",
                         "integrated.batch.png",
                         "integrated.lsi.clusters.png"]), category="Nuclear ATAC")
    params:
        outdir =lambda wildcards, output: dirname(output[0]),
        rscript= get_integrate_or_single_rscript, #join(ROOT_DIR, "R_scripts/annotations/integrateSignac.ipynb"), # The script defaults to the granja data
        sample_names = ",".join(samples.index),
        gff = gff
    shell: "papermill -p outdir {params.outdir} -p sample_names {params.sample_names} -p gff_id {params.gff} {params.rscript} {output[0]}"

#
# rule se_meta:
#     """Save seurat meta-data as csv"""
#     input:
#         se_f = "{outdir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds",
#     output:
#         se_meta = "{outdir}/data/annotation/gff_{gff}/mergedSamples/se_cells_meta.tsv",
#         note = "{outdir}/data/annotation/gff_{gff}/mergedSamples/se_cells_meta.ipynb",
#     params:
#         rscript = join(ROOT_DIR, "workflow/notebooks/lineage_clones/clusters_cells_meta.ipynb"),
#     shell: "papermill -p se_f {input.se_f} {params.rscript} {output.note}"
#
#
# def get_cluster_labels():
#     return config.get("umap_clusters_f", "FALSE")
#
# rule add_cluster_labels:
#     """Prepare clone-by-cluster counts for umap and hypergeometric test"""
#     input:
#         se_f = "{outdir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds",
#     output:
#         se_meta = "{outdir}/data/annotation/gff_{gff}/mergedSamples/se_cells_meta_labels.tsv",
#         note = "{outdir}/data/annotation/gff_{gff}/mergedSamples/add_cluster_labels.ipynb",
#     params:
#         labels = get_cluster_labels(),
#         rscript = join(ROOT_DIR, "workflow/notebooks/lineage_clones/add_cluster_labels.ipynb"),
#     shell: "papermill -p se_f {input.se_f} -p cluster_labels_f {params.labels} {params.rscript} {output.note}"
#

#######################################################################

def get_comps():
    if "comparisons" in config:
        return config["comparisons"]
    else:
        return 'NULL'


rule runDE:
    input:
        integrated = "{outdir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds",
    output:
        # report(expand(dir("{{outdir}}/data/annotationDE/{d}"),
        #               d=["clusters", "conditions_clusters", "conditions_conserved"])),
        "{outdir}/data/annotation/gff_{gff}/mergedSamples/DE/DE.ipynb",
        report(expand("{{outdir}}/data/annotation/gff_{{gff}}/mergedSamples/DE/{f}",
                f=["seuratImmuneDotPlot.png", "linImmuneDotPlot.png"]), category="Nuclear ATAC")
    params:
       # indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        rscript= join(ROOT_DIR, "R_scripts/annotations/samplesMergedSignac_DE.ipynb"), # The script defaults to the granja data
        sample_names = ",".join(samples.index),
        comps_f = get_comps()
        #samples = ",".join(samples["cellr_ID"])
    shell: "papermill -p integrated_f {input} -p outdir {params.outdir} -p sample_names {params.sample_names} -p comps_f {params.comps_f:q} {params.rscript} {output[0]}"


# rule runDE_rm_clust:
#     params:
#         clust_rm = config["clust_rm"]
#     output:
#         "{outdir}/data/annotation/gff_{gff}/mergedSamples/DE/DE.ipynb",
#     shell: "papermill -p clust_rm {params.clust_rm} -p integrated_f {input} -p outdir {params.outdir} -p sample_names {params.sample_names} -p comps_f {params.comps_f:q} {params.rscript} {output[0]}"


rule runDE_TF:
    input:
        integrated = "{outdir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds",
    output:
        note="{outdir}/data/annotation/gff_{gff}/mergedSamples/DE_TF/DE_TF.ipynb",
    params:
        outdir = lambda wildcards, output: dirname(output[0]),
        rscript= join(ROOT_DIR, "R_scripts/annotations/DE_TF.ipynb"), # The script defaults to the granja data
        sample_names = ",".join(samples.index),
        comps_f = get_comps(),
        genome= params["genome_path"][config["genome"]]["ref_fa"]
    shell: "papermill -p integrated_f {input} -p outdir {params.outdir} -p sample_names {params.sample_names} -p comps_f {params.comps_f:q} -p genome {params.genome} {params.rscript} {output[0]}"


rule runGSEA:
    input:
        DE_out_path = "{outdir}/data/annotation/gff_{gff}/mergedSamples/DE/DE.ipynb",
    output:
        note="{outdir}/data/annotation/gff_{gff}/mergedSamples/DE/GSEA/clusters/GSEA.ipynb",
    params:
        input = lambda wildcards, input: join(dirname(input[0]), "clusters"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "R_scripts/annotations/runGSEA.ipynb"),
        gsea_dir = join(ROOT_DIR, "software/Bioinformatics_Tools/"), #/ is necessary
    conda:
        "../envs/gsea.yml" #environment with clusterprofiler
    shell: "papermill -p DE.out.path {params.input} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"


# rule overlay_cells_meta:
#     input:
#         "{outdir}/data/annotationallSamples.integrated.rds",
#     output:
#         "{outdir}/data/annotationdonors.png",
#         "{outdir}/data/annotationclones.png"
#     params:
#         indir = lambda wildcards, input: dirname(input[0]),
#         outdir = lambda wildcards, output: dirname(output[0]),
#         cells_meta = config["cells_meta"],
#         rscript= join(ROOT_DIR, "R_scripts/annotations/lineageNuclear.ipynb"), # The script defaults to the granja data
#
#     shell: "papermill -p cellr_in {params.indir} -p outdir {params.outdir} -p cells_meta {params.cells_meta} {params.rscript} {output[0]}"
#
