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


rule all:
    input:
        expand("{outdir}/annotation/{prefix}/mergedSamples/allSamples.integrated.ipynb",
               outdir=config['outdir'], prefix=config["prefix"]),

        # expand("{outdir}/annotation/data/{extrnl}/{extrnl}.fragments.sort.tsv.gz",
        #        extrnl=config['annotations']['name'], outdir=config['outdir'],
        #        ),
        # expand("{outdir}/annotation/{prefix}/{sample}/lareau/{proj}/{sample}.clusters.csv",
        #        outdir=config['outdir'], sample=samples.index,
        #        proj=cfg_anno["lareau"]["proj_ref"],
        #        prefix=config["prefix"]),
         expand("{outdir}/annotation/{prefix}/mergedSamples/DE/seuratImmuneDotPlot.png",
                 outdir=config['outdir'],prefix=config["prefix"]),
         expand("{outdir}/annotation/{prefix}/mergedSamples/DE/GSEA/GSEA.ipynb",
                 outdir=config['outdir'],prefix=config["prefix"])

        # expand("{outdir}/annotation/{prefix}/{sample}/anchors/{assay}/{sample}.clusters.csv", sample=samples.index,
        #         assay=cfg_anno["anchors"]["assay"], outdir=config['outdir'],
        #         prefix=config["prefix"])

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
        "{outdir}/annotation/{prefix}/{sample}/{sample}.merged.ipynb",
        "{outdir}/annotation/{prefix}/{sample}/{sample}.merged.rds",
        report("{outdir}/annotation/{prefix}/{sample}/{sample}.merged.lsi.Batchlabels.png")
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


#conda-env-mito_trace-r
rule annotation_lareau:
    input: "{outdir}/annotation/{prefix}/{sample}/{sample}.merged.rds"
    output:
        report("{outdir}/annotation/{prefix}/{sample}/lareau/{proj}/{sample}.clusters.csv"), # "{outdir}/annotation/{sample}/{sample}.clusters.txt"
        report("{outdir}/annotation/{prefix}/{sample}/lareau/{proj}/{sample}.project.sample.cellLabels.abundace.png"),
        note="{outdir}/annotation/{prefix}/{sample}/lareau/{proj}/{sample}.clusters.ipynb"
    params:
        outdir=lambda wildcards, output: dirname(output[0]),
        exp = lambda wildcards: wildcards.sample,
        nTop = 25000,
        proj = lambda wildcards: wildcards.proj, #config["method"],
        rscript= join(ROOT_DIR, "R_scripts/annotations/02_predict_lareau.ipynb")
    shell: "papermill -k ir -p exp {params.exp} -p nTop {params.nTop} -p SE_f {input[0]} -p outdir {params.outdir}  {params.rscript} {output.note}"



rule annotation_anchors:
    input: "{outdir}/annotation/{prefix}/{sample}/{sample}.merged.rds"
    output:
        report("{outdir}/annotation/{prefix}/{sample}/anchors/{assay}/{sample}.clusters.csv"), # "{outdir}/annotation/{sample}/{sample}.clusters.txt"
        report("{outdir}/annotation/{prefix}/{sample}/anchors/{assay}/{sample}.merged.anchors.sample.labels.png"),
        report("{outdir}/annotation/{prefix}/{sample}/anchors/{assay}/{sample}.merged.anchors.labels.png"),
        note="{outdir}/annotation/{prefix}/{sample}/anchors/{assay}/{sample}.clusters.ipynb"
    params:
        outdir=lambda wildcards, output: dirname(output[0]),
        exp = lambda wildcards: wildcards.sample,
        assay = "RNA",
        rscript= join(ROOT_DIR, "R_scripts/annotations/02_predict_anchors.ipynb")
    shell: "papermill -p exp {params.exp} -p SE_f {input[0]} -p outdir {params.outdir} {params.rscript} {output[1]}"


rule createMergedExpSignac:
    """ Creates a merged R SummarizedExperiment object with the external data and reference combined.
    The peaks are reduced to the overlapping sets, and if a fragments file is provided,
    Will re-count peaks after creating overlapping peak set (and changing the lengths if the overlaps extend)

    """
    output:
        "{outdir}/annotation/{prefix}/mergedSamples/allSamples.integrated.ipynb",
        "{outdir}/annotation/{prefix}/mergedSamples/allSamples.integrated.rds",
        report(expand("{{outdir}}/annotation/{{prefix}}/mergedSamples/{f}",
                      f=["integrated.merged.compare.png", "integrated.batch.png", "integrated.lsi.clusters.png"]))
        #report("{outdir}/annotation/{prefix}/{sample}/{sample}.merged.lsi.Batchlabels.png")
    params:
        indir = config["mtscATAC_OUTDIR"],
        outdir =lambda wildcards, output: dirname(output[0]),
        rscript= join(ROOT_DIR, "R_scripts/annotations/samplesCreateMergedSignac.ipynb"), # The script defaults to the granja data
        sample_names = ",".join(samples.index),
        samples = ",".join(samples["cellr_ID"])
        #workdir = os.getcwd(),
    shell: "papermill -p cellr_in {params.indir} -p outdir {params.outdir} -p samples {params.samples} -p sample_names {params.sample_names} {params.rscript} {output[0]}"


def get_comps():
    if "comparisons" in config:
        return config["comparisons"]
    else:
        return 'NULL'

rule runDE:
    input:
        integrated = "{outdir}/annotation/{prefix}/mergedSamples/allSamples.integrated.rds",
    output:
        # report(expand(dir("{{outdir}}/annotation/{{prefix}}/DE/{d}"),
        #               d=["clusters", "conditions_clusters", "conditions_conserved"])),
        "{outdir}/annotation/{prefix}/mergedSamples/DE/DE.ipynb",
        report(expand("{{outdir}}/annotation/{{prefix}}/mergedSamples/DE/{f}",
                f=["seuratImmuneDotPlot.png", "linImmuneDotPlot.png"]))
    params:
       # indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        rscript= join(ROOT_DIR, "R_scripts/annotations/samplesMergedSignac_DE.ipynb"), # The script defaults to the granja data
        sample_names = ",".join(samples.index),
        comps_f = get_comps()
        #samples = ",".join(samples["cellr_ID"])
    shell: "papermill -p integrated_f {input} -p outdir {params.outdir} -p sample_names {params.sample_names} -p comps_f {params.comps_f:q} {params.rscript} {output[0]}"


rule runGSEA:
    input:
        DE_out_path = "{outdir}/annotation/{prefix}/mergedSamples/DE/DE.ipynb",
    output:
        "{outdir}/annotation/{prefix}/mergedSamples/DE/GSEA/GSEA.ipynb",
    params:
        output = lambda wildcards, output: dirname(output[0]),
        rscript = join(ROOT_DIR, "R_scripts/annotations/runGSEA.ipynb")
    shell: "papermill -p DE.out.path {input} -p export.path {params.output} {params.rscript} {output[0]}"


# rule overlay_cells_meta:
#     input:
#         "{outdir}/annotation/{prefix}/allSamples.integrated.rds",
#     output:
#         "{outdir}/annotation/{prefix}/donors.png",
#         "{outdir}/annotation/{prefix}/clones.png"
#     params:
#         indir = lambda wildcards, input: dirname(input[0]),
#         outdir = lambda wildcards, output: dirname(output[0]),
#         cells_meta = config["cells_meta"],
#         rscript= join(ROOT_DIR, "R_scripts/annotations/lineageNuclear.ipynb"), # The script defaults to the granja data
#
#     shell: "papermill -p cellr_in {params.indir} -p outdir {params.outdir} -p cells_meta {params.cells_meta} {params.rscript} {output[0]}"
#
