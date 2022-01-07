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
#         # expand("{outdir}/annotation/data/{extrnl}/{extrnl}.fragments.sort.tsv.gz",
#         #        extrnl=config['annotations']['name'], outdir=config['outdir'],
#         #        ),
#         # expand("{outdir}/data/annotation{sample}/lareau/{proj}/{sample}.clusters.csv",
#         #        outdir=config['outdir'], sample=samples.index,
#         #        proj=cfg_anno["lareau"]["proj_ref"],
#         #        prefix=config["prefix"]),
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


#conda-env-mito_trace-r
rule annotation_lareau:
    input: "{outdir}/data/annotation{sample}/{sample}.merged.rds"
    output:
        report("{outdir}/data/annotation{sample}/lareau/{proj}/{sample}.clusters.csv"), # "{outdir}/annotation/{sample}/{sample}.clusters.txt"
        report("{outdir}/data/annotation{sample}/lareau/{proj}/{sample}.project.sample.cellLabels.abundace.png"),
        note="{outdir}/data/annotation{sample}/lareau/{proj}/{sample}.clusters.ipynb"
    params:
        outdir=lambda wildcards, output: dirname(output[0]),
        exp = lambda wildcards: wildcards.sample,
        nTop = 25000,
        proj = lambda wildcards: wildcards.proj, #config["method"],
        rscript= join(ROOT_DIR, "R_scripts/annotations/02_predict_lareau.ipynb")
    shell: "papermill -k ir -p exp {params.exp} -p nTop {params.nTop} -p SE_f {input[0]} -p outdir {params.outdir}  {params.rscript} {output.note}"


rule annotation_anchors:
    input: "{outdir}/data/annotation{sample}/{sample}.merged.rds"
    output:
        report("{outdir}/data/annotation{sample}/anchors/{assay}/{sample}.clusters.csv"), # "{outdir}/annotation/{sample}/{sample}.clusters.txt"
        report("{outdir}/data/annotation{sample}/anchors/{assay}/{sample}.merged.anchors.sample.labels.png"),
        report("{outdir}/data/annotation{sample}/anchors/{assay}/{sample}.merged.anchors.labels.png"),
        note="{outdir}/data/annotation{sample}/anchors/{assay}/{sample}.clusters.ipynb"
    params:
        outdir=lambda wildcards, output: dirname(output[0]),
        exp = lambda wildcards: wildcards.sample,
        assay = "RNA",
        rscript= join(ROOT_DIR, "R_scripts/annotations/02_predict_anchors.ipynb")
    shell: "papermill -p exp {params.exp} -p SE_f {input[0]} -p outdir {params.outdir} {params.rscript} {output[1]}"
