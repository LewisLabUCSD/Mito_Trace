from src.config import ROOT_DIR
import pandas as pd
from os.path import join, dirname
import os
import numpy as np
from icecream import ic


wildcard_constraints:
    variants = "simple|mgatkdonor",
    d = "%d"


samples = config["samples"] #pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)

clones_cfg = config["params"]["clones"]
nclonelist = clones_cfg['vireo']['params']['nclonelist']

########################################################################
## Variant Intermediate step for other workflows.
## Extract pileups for each donor from original filter matrices and run mgatk
########################################################################
# def get_coverage(wildcards):
#     return f"{config['cov_indir']['sample']}/{wildcards.sample}.coverage.txt"
def concat_knn(in_files, out_name):
    all = []
    cols = ["donor", "lineage"]
    for i in in_files:
        #print(i)
        all.append(pd.read_csv(i,sep='\t'))
        #print(all[-1].head())
        if "donor_index" in all[-1].columns.values:
            cols.append("donor_index")
        if "lineage_index" in all[-1].columns.values:
            cols.append("lineage_index")
    all = pd.concat(all,ignore_index=True).sort_values(cols)
    if 'level_0' in all.columns:
        all = all.drop('level_0',axis=1)
    all.to_csv(out_name,sep='\t',index=False)


def get_knn_script(cfg):
    if cfg["variants_type"] == "init":
        join(ROOT_DIR, "R_scripts", "knn_clones_init.ipynb"),
    elif cfg["variants_type"] == "mgatkdonor":
        return join(ROOT_DIR, "R_scripts", "knn_clones.ipynb"),

rule knn:
    input:
        af ="{outdir}/donor{d}/af.tsv"
    output:
        cells= "{outdir}/knn/kparam_{kparam}/donors/donor{d}/cells_meta.tsv",
        fig = "{outdir}/knn/kparam_{kparam}/donors/donor{d}/donor{d}.variants.labels.png"
    params:
        indir = lambda wildcards, input: dirname(input.af),
        name = lambda wildcards: f"donor{wildcards.d}",
        outdir = lambda wildcards, output: dirname(output[0]),
        note = get_knn_script(config) #join(ROOT_DIR, "R_scripts", "knn_clones_init.ipynb"),
    shell:
        "papermill -k ir -p indir {params.indir} -p name {params.name} -p donor {wildcards.d} -p outdir {params.outdir} -p kparam {wildcards.kparam} {params.note} {params.outdir}/output.ipynb"


rule knn_concat:
    input:
        cells= expand("{{outdir}}/knn/kparam_{{kparam}}/donors/donor{d}/cells_meta.tsv",
                       d=np.arange(config["N_DONORS"])),
    output:
        cells_all="{outdir}/knn/kparam_{kparam}/cells_meta.tsv"
    run: concat_knn(input, output[0])


# Sort of like an all rule to bring downstream results together.
rule complete_lineage:
    input:
        "{outdir}/clones/variants_{variants}/{method}/temp/.tmp",
    output: "{outdir}/clones/variants_{variants}/{method}/.completed"
    shell: "touch {output}"


