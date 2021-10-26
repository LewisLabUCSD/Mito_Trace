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
        expand("{outdir}/annotation_clones/data/{prefix}/SE.rds",
               outdir=config['outdir'], prefix=config["prefix"]),
        expand("{outdir}/annotation_clones/data/{prefix}/DE/donor{n}/DE.ipynb",
                outdir=config['outdir'], prefix=config["prefix"],
                n=config['N_DONORS']),
        expand("{outdir}/annotation_clones/data/{prefix}/DE_TF/donor{n}/DE_TF.ipynb",
                outdir=config['outdir'], prefix=config["prefix"],
                n=config['N_DONORS']),

        expand("{outdir}/annotation_clones/data/{prefix}/DE/donor{n}/GSEA/clusters/GSEA.ipynb",
                outdir=config['outdir'], prefix=config["prefix"],
                n=config['N_DONORS']),

        expand("{outdir}/annotation_clones/data/{prefix}/allDonors/DE/DE.ipynb",
                outdir=config['outdir'], prefix=config["prefix"],
                ),


def get_comps():
    if "comparisons" in config:
        return config["comparisons"]
    else:
        return 'NULL'


rule addClones:
    input:
        noc = "{outdir}/annotation_clones/data/{prefix}/mergedSamples/integrated.rds",
        clones = "{indir}/cells_meta.tsv",#"{outdir}/pipeline/published/{prefix}/data/clones/clones.txt",
    output:
        se_f = "{outdir}/annotation_clones/data/{prefix}/SE.rds",
        note = "{outdir}/annotation_clones/data/{prefix}/addClones.ipynb"
    params:
        outdir = lambda wildcards, output: dirname(output[0]),
        rscript= join(ROOT_DIR, "R_scripts/annotations/clones/addClones_01.vCurrent.ipynb"), # The script defaults to the granja data
    shell: "papermill -p cells_meta_f {input.clones} -p se_f {input.noc} -p outdir {params.outdir} {params.rscript} {output.note}"


rule runDE_enrich:
    input:
        se_f = "{outdir}/annotation_clones/data/{prefix}/SE.rds",
        enrich_f = "{enrich_indir}/cells_meta.tsv",
    output:
        "{outdir}/annotation_clones/data/{prefix}/DE/donor{d}/DE.ipynb"
    params:
        d = lambda wildcards: wildcards.d,
        outdir = lambda wildcards, output: dirname(output),
        rscript= join(ROOT_DIR, "R_scripts/annotations/clones/addEnrichment.vCurrent.ipynb")
    shell: "papermill -p enrich_f {input.enrich_f} -p se_f {input.se_f} -p outdir {params.outdir} {params.rscript} {output.note}"

# rule runDE:
#     input:
#         clones = "{outdir}/annotation_clones/data/{prefix}/SE.rds"
#     output:
#         "{outdir}/annotation_clones/data/{prefix}/DE/donor{n}/DE.ipynb",
#     params:
#         outdir = lambda wildcards, output: dirname(output[0]),
#         rscript= join(ROOT_DIR, "R_scripts/annotations/clones/DE_clones.vCurrent.ipynb"), # The script defaults to the granja data
#         sample_names = ",".join(samples.index),
#         comps_f = get_comps()
#         #samples = ",".join(samples["cellr_ID"])
#     shell: "papermill -p se {input} -p outdir {params.outdir} -p sample_names {params.sample_names} -p comps_f {params.comps_f:q} {params.rscript} {output[0]}"


rule runGSEA:
    input:
        DE_out_path = "{outdir}/annotation_clones/data/{prefix}/DE/donor{n}/DE.ipynb"
    output:
        "{outdir}/annotation_clones/data/{prefix}/DE/donor{n}/GSEA/clusters/GSEA.ipynb",
    params:
        input = lambda wildcards, input: join(dirname(input[0]), "clusters"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "R_scripts/annotations/clones/runGSEA_clones.ipynb")
    shell: "papermill -p DE.out.path {params.input} -p export.path {params.output} {params.rscript} {output[0]}"


# rule runDE_TF:
#     input:
#         integrated = "{outdir}/annotation_clones/data/{prefix}/{prefix}/mergedSamples/allSamples.integrated.rds",
#     output:
#         "{outdir}/annotation_clones/data/{prefix}/{prefix}/mergedSamples/DE_TF/DE_TF.ipynb",
#     params:
#         outdir = lambda wildcards, output: dirname(output[0]),
#         rscript= join(ROOT_DIR, "R_scripts/annotations/DE_TF.ipynb"), # The script defaults to the granja data
#         sample_names = ",".join(samples.index),
#         comps_f = get_comps()
#     shell: "papermill -p integrated_f {input} -p outdir {params.outdir} -p sample_names {params.sample_names} -p comps_f {params.comps_f:q} {params.rscript} {output[0]}"
#


# rule overlay_cells_meta:
#     input:
#         "{outdir}/annotation_clones/data/{prefix}/{prefix}/allSamples.integrated.rds",
#     output:
#         "{outdir}/annotation_clones/data/{prefix}/{prefix}/donors.png",
#         "{outdir}/annotation_clones/data/{prefix}/{prefix}/clones.png"
#     params:
#         indir = lambda wildcards, input: dirname(input[0]),
#         outdir = lambda wildcards, output: dirname(output[0]),
#         cells_meta = config["cells_meta"],
#         rscript= join(ROOT_DIR, "R_scripts/annotations/lineageNuclear.ipynb"), # The script defaults to the granja data
#
#     shell: "papermill -p cellr_in {params.indir} -p outdir {params.outdir} -p cells_meta {params.cells_meta} {params.rscript} {output[0]}"
#
