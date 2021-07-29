#configfile: "parameters/Zhu_Single_Cell_10X_genomics_Human2_002.yaml"
from os.path import join, dirname
from src.utils.data_io import sparse_to_cellranger_fragments
import pandas as pd
from src.config import ROOT_DIR

#configfile: "/data2/mito_lineage/parameters/DUPI_april08_2021/mttrace_mtasnucl.yaml"

cfg_anno = config['annotations']
extrnl = cfg_anno["name"]
samples = pd.read_table(join(ROOT_DIR, config["samples"]), dtype=str,sep=',').set_index(["sample_name"], drop=False)
res = config["results"]


rule all:
    input:
        expand("output/data/{extrnl}/{extrnl}.fragments.tsv", extrnl=config['annotations']['name']),
        expand("output/{results}/{sample}/{sample}.clusters.txt", results=res, sample=samples['sample_name'].values)
        #f"data/processed/external/{extrnl}/{extrnl}.fragments.tsv"#, out=config['annotations']['name'])


def get_cellr_input(wildcards):
    # peaks_f = f"{indir}/peaks.bed"
    # singlecell_f = f"{indir}/single.bed"
    # mtx_f = ""
    out_d = cfg_anno['datadir']
    return join(out_d, cfg_anno['mtx_f']), join(out_d, cfg_anno['peaks_f']), join(out_d, cfg_anno['cells_f'])

rule peaks_to_fragments:
    input: get_cellr_input
    output: "output/data/{extrnl}/{extrnl}.fragments.tsv"
    run: sparse_to_cellranger_fragments(input[0], input[1], input[2], out_f=output[0], n_cpus=24)

# rule project_annotation:
#     input: "/data2/mito_lineage/data/processed/mttrace/{prefix}/MTblacklist/{sample}/MT/cellr_True/{sample}_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/"
#     output:  "/data2/mito_lineage/Analysis/annotation/output/data/{sample}"
#     params:
#         name= lambda wildcards: wildcards.sample,
#         script= join("src", "R_scripts", "annotations", "annotations.ipynb")
#     shell: "papermill -p name {params.name} -p params.indir {params.indir} -p outdir {params.outdir} -in_mgatk {params.in_mgatk} {output[0]}"


def get_sample_barcodes(wildcards):
    print(samples)
    print(wildcards)
    return samples.loc[wildcards.sample, "barcode_f"]


rule createSE:
    input: get_sample_barcodes
    output: "output/{results}/{sample}/{sample}.merged.rds"
    params:
        outdir=lambda wildcards, output: dirname(output[0]),
        indir=lambda wildcards, input: dirname(input[0]),
        rscript= join(ROOT_DIR, "R_scripts/annotations/create_SE.ipynb"), # The script defaults to the granja data
        exp = lambda wildcards: wildcards.sample,
    shell: "Rscript R_scripts/annotations/create_SE.ipynb -p singlecell_sumstats_dir {params.indir} -p exp {params.exp} -p outdir {params.outdir}"


rule annotation_peaks:
    input: "output/{results}/{sample}/{sample}.merged.rds"
    output: "output/{results}/{sample}/{sample}.clusters.txt" # "output/{results}/{sample}/{sample}.clusters.txt"
    params:
        outdir=lambda wildcards, output: dirname(output[0]),
        exp = lambda wildcards: wildcards.sample,
        nTop = 25000,
        rscript= join(ROOT_DIR, "R_scripts/annotations/orig_01_CD34_projection.ipynb")
    shell: "Rscript {params.rscript} -p exp {params.exp} -p nTop {params.nTop} -p outdir {params.outdir}"
