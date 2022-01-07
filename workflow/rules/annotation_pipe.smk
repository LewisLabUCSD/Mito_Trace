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


# rule all:
#     input:
#         expand("{outdir}/data/annotation/mergedSamples/allSamples.integrated.ipynb",
#                outdir=config['outdir'], prefix=config["prefix"]),
#
#         # expand("{outdir}/annotation/data/{extrnl}/{extrnl}.fragments.sort.tsv.gz",
#         #        extrnl=config['annotations']['name'], outdir=config['outdir'],
#         #        ),
#         # expand("{outdir}/data/annotation{sample}/lareau/{proj}/{sample}.clusters.csv",
#         #        outdir=config['outdir'], sample=samples.index,
#         #        proj=cfg_anno["lareau"]["proj_ref"],
#         #        prefix=config["prefix"]),
#          expand("{outdir}/data/annotation/mergedSamples/DE/seuratImmuneDotPlot.png",
#                  outdir=config['outdir'],prefix=config["prefix"]),
#
#         expand("{outdir}/data/annotation/mergedSamples/DE_TF/DE_TF.ipynb",
#                  outdir=config['outdir'],prefix=config["prefix"]),
#          expand("{outdir}/data/annotation/mergedSamples/DE/GSEA/clusters/GSEA.ipynb",
#                  outdir=config['outdir'],prefix=config["prefix"])
#
#         # expand("{outdir}/data/annotation{sample}/anchors/{assay}/{sample}.clusters.csv", sample=samples.index,
#         #         assay=cfg_anno["anchors"]["assay"], outdir=config['outdir'],
#         #         prefix=config["prefix"])

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


rule createMergedSignac:
    """ Creates a merged R SummarizedExperiment object with the external data and reference combined.
    The peaks are reduced to the overlapping sets, and if a fragments file is provided,
    Will re-count peaks after creating overlapping peak set (and changing the lengths if the overlaps extend)

    """
    input:
        sample_indir=get_cellr_dir
        #"{outdir}/annotation/data/{extrnl}/{extrnl}.fragments.sort.tsv.gz"

    output:
        "{outdir}/data/annotation{sample}/{sample}.merged.ipynb",
        "{outdir}/data/annotation{sample}/{sample}.merged.rds",
        report("{outdir}/data/annotation{sample}/{sample}.merged.lsi.Batchlabels.png")
    params:
        outdir=lambda wildcards, output: dirname(output[0]),
        exp = lambda wildcards: wildcards.sample,
        #cells_subset=None, To add. Text file of cell labels.
        # # External
        external_frag_file = get_extrnl_frags, #lambda wc: f"{wc.outdir}/annotation/data/{wc.outdir}/{wc.outdir}.fragments.sort.tsv.gz",
        external_prefix = cfg_anno["ID"],
        rscript= join(ROOT_DIR, "R_scripts/annotations/01_createMergedSignac.ipynb"), # The script defaults to the granja data
        #workdir = os.getcwd(),
    shell: "papermill  -p sample_indir {input} -p exp {params.exp} -p outdir {params.outdir} -p external_prefix {params.external_prefix} -p external_frag_file {params.external_frag_file} {params.rscript} {output[0]}"


#####################
## V02
#####################
rule createExpSignac:
    """ Creates a merged R SummarizedExperiment object with the external data and reference combined.
    The peaks are reduced to the overlapping sets, and if a fragments file is provided,
    Will re-count peaks after creating overlapping peak set (and changing the lengths if the overlaps extend)
    """
    output:
        "{outdir}/data/annotation/mergedSamples/createExpSignac.ipynb",
        "{outdir}/data/annotation/mergedSamples/allSamples.rds",
        report("{outdir}/data/annotation/mergedSamples/QC_01.png", category="Nuclear ATAC",
                subcategory="QC"),
        "{outdir}/data/annotation/mergedSamples/QC_01.pdf"
    params:
        indir = config["mtscATAC_OUTDIR"],
        outdir =lambda wildcards, output: dirname(output[0]),
        rscript= join(ROOT_DIR, "R_scripts/annotations/samplesCreateSignac.ipynb"),
        sample_names = ",".join(samples.index),
        samples = ",".join(samples["cellr_ID"].values)
        #workdir = os.getcwd(),
    shell: "papermill -p cellr_in {params.indir} -p outdir {params.outdir} -p samples {params.samples} -p sample_names {params.sample_names} {params.rscript} {output[0]}"


rule integrateSignac:
    input:
        a="{outdir}/data/annotation/mergedSamples/allSamples.rds",
    output:
        a="{outdir}/data/annotation/mergedSamples/allSamples.integrated.ipynb",
        b="{outdir}/data/annotation/mergedSamples/allSamples.integrated.rds",
        c=report(expand("{{outdir}}/data/annotation/mergedSamples/{f}",
                      f=["QC_02.png",
                         "integrated.merged.compare.png",
                         "integrated.batch.png",
                         "integrated.lsi.clusters.png"]), category="Nuclear ATAC")
    params:
        outdir =lambda wildcards, output: dirname(output[0]),
        rscript= join(ROOT_DIR, "R_scripts/annotations/integrateSignac.ipynb"), # The script defaults to the granja data
        sample_names = ",".join(samples.index),
    shell: "papermill -p outdir {params.outdir} -p sample_names {params.sample_names} {params.rscript} {output[0]}"
#######################################################################

def get_comps():
    if "comparisons" in config:
        return config["comparisons"]
    else:
        return 'NULL'


rule runDE:
    input:
        integrated = "{outdir}/data/annotation/mergedSamples/allSamples.integrated.rds",
    output:
        # report(expand(dir("{{outdir}}/data/annotationDE/{d}"),
        #               d=["clusters", "conditions_clusters", "conditions_conserved"])),
        "{outdir}/data/annotation/mergedSamples/DE/DE.ipynb",
        report(expand("{{outdir}}/data/annotation/mergedSamples/DE/{f}",
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
#         "{outdir}/data/annotation/mergedSamples/DE/DE.ipynb",
#     shell: "papermill -p clust_rm {params.clust_rm} -p integrated_f {input} -p outdir {params.outdir} -p sample_names {params.sample_names} -p comps_f {params.comps_f:q} {params.rscript} {output[0]}"


rule runDE_TF:
    input:
        integrated = "{outdir}/data/annotation/mergedSamples/allSamples.integrated.rds",
    output:
        note="{outdir}/data/annotation/mergedSamples/DE_TF/DE_TF.ipynb",
    params:
        outdir = lambda wildcards, output: dirname(output[0]),
        rscript= join(ROOT_DIR, "R_scripts/annotations/DE_TF.ipynb"), # The script defaults to the granja data
        sample_names = ",".join(samples.index),
        comps_f = get_comps(),
        genome= params["genome_path"][config["genome"]]["ref_fa"]
    shell: "papermill -p integrated_f {input} -p outdir {params.outdir} -p sample_names {params.sample_names} -p comps_f {params.comps_f:q} -p genome {params.genome} {params.rscript} {output[0]}"


rule runGSEA:
    input:
        DE_out_path = "{outdir}/data/annotation/mergedSamples/DE/DE.ipynb",
    output:
        note="{outdir}/data/annotation/mergedSamples/DE/GSEA/clusters/GSEA.ipynb",
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
