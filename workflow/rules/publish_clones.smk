
#report: "report/workflow.rst"
from os.path import dirname
from os.path import join
#from src.config import ROOT_DIR
#workdir: ROOT_DIR
wildcard_constraints:
    cellr='True|False',
    kparam='[0-9]+',  #'%d', #[0-9]+',
    variants = "simple|mgatkdonor|init|simpleUnion",
    d = "[0-9]+"

#configfile:  "parameters/pipeline/cosmo_server/jan21_2021.yaml"
from src.config import ROOT_DIR
from src.utils.parse_config import read_config_file
import os
import numpy as np
import pandas as pd
from snakemake.utils import min_version
from icecream import ic
min_version("6.0")
print('config', config)

########################################################################
# Setup parameters and outdir
########################################################################
res = join(config["outdir"], "pipeline", config["prefix"])

########################
# Setup parameters for each analysis
########################
params = read_config_file(config["config"])
version = params["version"]

samples = pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)
anno_res = join(config["outdir"], "annotation", "data", config["prefix"])

params_mt = params["mtpreproc"]
params_ft = params["filters"]
ft = params_ft["params"]
params_mult = params["multiplex"]
params_mgatk = params["mgatk"]
params_clones = params["clones"]
params_annclo = params["annotation_clones"]["params"]
params_clch = params["de_clones_change"]["params"]
params_clch["samples"] = samples


ref_fa = params["genome_path"][config["genome"]]["ref_fa"]
mt_ref_fa = params["genome_path"][config["genome"]]["mt_ref_fa"]



########################
# Define variables here:
########################
out_dir=res
gff = params["genome_path"][config["genome"]]["gff"]
cellrbc = params_mt["params"]["cellrbc"]
num_read = params_mt["params"]["numreadsfilter"]
maxBP = params_mt["maxBP"]
mincells=ft['mincells']
minreads=ft['minreads']
topN=ft["topN"]
hetthresh=ft['hetthresh']
minhetcells=ft['minhetcells']
hetcountthresh=ft['hetcountthresh']
bqthresh=ft['bqthresh']
kparam=params_clones["knn"]["params"]["resolution"]
minPct=params_annclo["min_pct"]
gsea_pval=params_annclo["gsea_pval"]
min_cells=params_annclo["min_cells"]
variants=params_clones["variants"]
stat_col = params_annclo["stat_col"]
padjmethod = params_annclo["padjmethod"]
filt = params_clch["filt"]
shuffle = params_clch["shuffle"]
use_p_adjust = params_clch["use_p_adjust"]
clch_pthresh = params_clch["pthresh"]
min_cell=params_clch["min_cell"]
min_cell_both=params_clch["min_cell_both"]

#hyperg=expand("{{outdir}}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/hypergeom.csv",
hyperMinCl=params_annclo["hyperMinCl"]
bothConds=params_annclo["bothConds"]
pthresh=params_annclo["p_thresh"]



def create_clone_report(inputs, output):
    return

###################################
## Rules and workflow
###################################
###################################
rule all:
    input:
        # expand("{out_dir}/reports/summary_clones/v{id}/clones.pdf",
        #     out_dir=res, id=params["version"]),
        expand("{out_dir}/reports/summary_clones/v{id}/clones.ipynb", out_dir=res, id=version)

rule clone_report:
    input:
        anno_dir = f"{out_dir}/data/annotation/gff_{gff}/mergedSamples/",
        clone_dir = f"{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/",
        annClo_dir = f"{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/",
        aggreg_dir = f"{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{clch_pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/merge_enrich_and_lineage",
        hypergeom_dir = f"{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}",
    output:
       # pdf = "{out_dir}/reports/summary_clones/v{id}/clones.pdf",
        note = "{out_dir}/reports/summary_clones/v{id}/clones.ipynb",
    params:
        rscript = "workflow/notebooks/summarize_clones.ipynb",
        params_f = config["config"],
        outdir = lambda wildcards, output: dirname(output.note)
    shell:
        "papermill -p outdir {params.outdir} -p anno_dir {input.anno_dir} -p clone_dir {input.clone_dir} -p aggreg_dir {input.aggreg_dir} -p annClo_dir {input.annClo_dir} -p hypergeom {input.hypergeom_dir} -p params_f {params.params_f} {params.rscript} {output.note}"
        #create_clone_report(input, output)

# rule rclone:
#     input:
#         "{input}"
#     output:
#         "{input}/output/rsync_copy.ok"
#     params:
#         rsync_loc = config["rsync_location"]
#     shell: "rcopy {input} drive:{rsync_loc} > rsync_copy.ok && mv rsync_copy.ok {output}/"
