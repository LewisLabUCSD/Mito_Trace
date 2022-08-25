import pandas as pd
from snakemake.utils import min_version
from icecream import ic
min_version("6.0")
print('config', config)
from src.config import ROOT_DIR
from os.path import join, dirname
import numpy as np
import pickle


#dendro_d = f"dendro_dt_{config['clones']['params']['dendro_thresh']}"

dendro_d = config['clones']['params']['dendro_thresh']

# inputOnly and noInput:
# See if Input there, then run inputOnly, and if other conditions also there, run noInput
sample_ids = config["samples"].index
print('sample_ids', sample_ids)
condition = []
if "Input" in sample_ids:
    condition.append("inputOnly")
    if len(sample_ids) > 1:
        condition.append("noInput")
elif len(sample_ids) > 0:
    condition.append("noInput")

if len(condition) == 0:
    raise ValueError("Samples not set up properly")
print('condition types for clonal shift', condition)


rule clonalshift_clones:
    input:
        se_meta = expand("{{outdir}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/gff_{gff}/annotation_clones/se_cells_meta_labels.tsv",
                         gff=config["gff"])
    output:
        note = "{outdir}/clonal_shifts/variants_{variants}/clones/results/{condition}/knn/kparam_{kparam}/output.ipynb"
    params:
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/hypergeometric_clones.ipynb"), #get_script,
        outdir = lambda wildcards, output: dirname(output.note),
        N_DONORS = config["N_DONORS"],
    shell: "papermill -p se_cells_meta_f {input.se_meta} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} -p condition {wildcards.condition} {params.script} {output}"



rule clonalshift_dendro_bc:
    input:
        se_meta = expand("{{outdir}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/gff_{gff}/annotation_clones/se_cells_meta_labels.tsv",
                         gff=config["gff"]),
        barcodes_dir = expand("{{outdir}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/barcodes/btwnClones_dendro_dt_{dt}/donor{d}.clones_dendro.csv",
                              dt=dendro_d, d=np.arange(config["N_DONORS"]))
        #barcodes_dir = "{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/btwnClones_dendro_dt_{dt}",
    output:
        note = "{outdir}/clonal_shifts/variants_{variants}/dendro_bc/results/{condition}/knn/kparam_{kparam}/output.ipynb"
    params:
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/hypergeometric_dendro.ipynb"), #get_script,
        outdir = lambda wildcards, output: dirname(output.note),
        N_DONORS = config["N_DONORS"],
        barcodes_dir = lambda wildcards, input: dirname(input.barcodes_dir[0]),
    shell: "papermill -p se_cells_meta_f {input.se_meta} -p outdir {params.outdir} -p barcodes_dir {params.barcodes_dir} -p N_DONORS {params.N_DONORS} -p condition {wildcards.condition} {params.script} {output}"

#
# rule compare_input_and_culture_dendro_bc:
#     input:
#         noInput = "{outdir}/clonal_shifts/variants_{variants}/dendro_bc/results/noInput/knn/kparam_{kparam}/output.ipynb",
#         input = "{outdir}/clonal_shifts/variants_{variants}/dendro_bc/results/inputOnly/knn/kparam_{kparam}/output.ipynb"
#     output:
#         note = "{outdir}/clonal_shifts/variants_{variants}/donors/donor{d}/dendro_bc/knn_kparam_{kparam}/output.ipynb"
#     params:
#         noInput_indir= lambda wildcards, input: dirname(input.noInput),
#         input_indir = lambda wildcards, input: dirname(input.input),
#         outdir = lambda wildcards, output: dirname(output.note),
#         donor = lambda wildcards: wildcards.d,
#         script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/combine_conditions_hypergeometric.ipynb"), #get_script,
#     shell: "papermill -p noIn_indir {params.noInput_indir} -p input_indir {params.input_indir} -p donor {params.donor} -p clone_col den_clust -p outdir {params.outdir} {params.script} {output.note}"
#


rule clonalshift_mt_as_clones:
    input:
        clones = "{outdir}/mt_as_clones/variants_{variants}/bestparams/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/best_params_save.ipynb",
        se_meta = expand("{{outdir}}/anno_multiplex/gff_{gff}/se_cells_meta_labels.tsv", gff=config["gff"])
    output:
        "{outdir}/clonal_shifts/variants_{variants}/mt_as_clones/results/{condition}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
    params:
        indir = lambda wildcards, input: dirname(input.clones),
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/hypergeometric_mt_as_clones.ipynb"), #get_script,
        outdir = lambda wildcards, output: dirname(output[0]),
        N_DONORS = config["N_DONORS"],
    shell: "papermill -p indir {params.indir} -p se_cells_meta_f {input.se_meta} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} -p condition {wildcards.condition} {params.script} {output}"


# rule compare_input_and_culture_mt_as_clones:
#     input:
#         noInput = "{outdir}/clonal_shifts/variants_{variants}/mt_as_clones/results/noInput/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
#         input = "{outdir}/clonal_shifts/variants_{variants}/mt_as_clones/results/inputOnly/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
#     output:
#         note = "{outdir}/clonal_shifts/variants_{variants}/donors/donor{d}/mt_as_clones/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
#     params:
#         noInput_indir= lambda wildcards, input: dirname(input.noInput),
#         input_indir = lambda wildcards, input: dirname(input.input),
#         outdir = lambda wildcards, output: dirname(output.note),
#         donor = lambda wildcards: wildcards.d,
#         script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/combine_conditions_hypergeometric.ipynb"), #get_script,
#     shell: "papermill -p noIn_indir {params.noInput_indir} -p input_indir {params.input_indir} -p donor {params.donor} -p clone_col Variants -p outdir {params.outdir} {params.script} {output.note}"
#

rule clonalshift_mt_as_clones_dendro:
    input:
        clones = "{outdir}/mt_as_clones/variants_{variants}/dendro/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/dendro_mt_clones.ipynb",
        se_meta = expand("{{outdir}}/anno_multiplex/gff_{gff}/se_cells_meta_labels.tsv", gff=config["gff"]) #"{outdir}/annotation_clones/se_cells_meta_labels.tsv"
    output:
        "{outdir}/clonal_shifts/variants_{variants}/mt_as_clones_dendro/results/{condition}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
    params:
        indir = lambda wildcards, input: dirname(input.clones),
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/hypergeometric_mt_as_clones_dendro.ipynb"), #get_script,
        outdir = lambda wildcards, output: dirname(output[0]),
        N_DONORS = config["N_DONORS"]
    shell: "papermill -p indir {params.indir} -p se_cells_meta_f {input.se_meta} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} -p condition {wildcards.condition} {params.script} {output}"


# rule compare_input_and_culture_mt_as_clones_dendro:
#     input:
#         noInput = "{outdir}/clonal_shifts/variants_{variants}/mt_as_clones_dendro/results/noInput/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
#         input = "{outdir}/clonal_shifts/variants_{variants}/mt_as_clones_dendro/results/inputOnly/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
#     output:
#         note = "{outdir}/clonal_shifts/variants_{variants}/donors/donor{d}/mt_as_clones_dendro/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
#     params:
#         noInput_indir= lambda wildcards, input: dirname(input.noInput),
#         input_indir = lambda wildcards, input: dirname(input.input),
#         outdir = lambda wildcards, output: dirname(output.note),
#         donor = lambda wildcards: wildcards.d,
#         script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/combine_conditions_hypergeometric.ipynb"), #get_script,
#     shell: "papermill -p noIn_indir {params.noInput_indir} -p input_indir {params.input_indir} -p donor {params.donor} -p clone_col den_clust -p outdir {params.outdir} {params.script} {output.note}"


name_map = {"clones":"name", "dendro_bc": "den_clust", "mt_as_clones": "Variants", "mt_as_clones_dendro": "den_clust"}

rule compare_input_and_culture_cl:
    input:
        noInput = "{outdir}/clonal_shifts/variants_{variants}/{cloneShift_method}/results/noInput/knn/kparam_{kparam}/output.ipynb",
        input = "{outdir}/clonal_shifts/variants_{variants}/{cloneShift_method}/results/inputOnly/knn/kparam_{kparam}/output.ipynb"
    output:
        note = "{outdir}/clonal_shifts/variants_{variants}/donors/donor{d}/{cloneShift_method}/knn_kparam_{kparam}/output_compare.ipynb"
    params:
        noInput_indir= lambda wildcards, input: dirname(input.noInput),
        input_indir = lambda wildcards, input: dirname(input.input),
        outdir = lambda wildcards, output: dirname(output.note),
        donor = lambda wildcards: wildcards.d,
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/combine_conditions_hypergeometric.ipynb"), #get_script,
        clone_col = lambda wildcards: name_map[wildcards.cloneShift_method]
    shell: "papermill -p noIn_indir {params.noInput_indir} -p input_indir {params.input_indir} -p donor {params.donor} -p clone_col {params.clone_col} -p outdir {params.outdir} {params.script} {output.note}"


rule single_input_and_culture_cl:
    input:
        cond = "{outdir}/clonal_shifts/variants_{variants}/{cloneShift_method}/results/{condition}/knn/kparam_{kparam}/output.ipynb",
    output:
        note = "{outdir}/clonal_shifts/variants_{variants}/donors/donor{d}/{cloneShift_method}/knn_kparam_{kparam}/output_single_{condition}_compare.ipynb"
    params:
        indir= lambda wildcards, input: dirname(input.cond),
        outdir = lambda wildcards, output: dirname(output.note),
        donor = lambda wildcards: wildcards.d,
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/individual_conditions_hypergeometric.ipynb"), #get_script,
        cond = condition[0],# inputOnly or noInput
        clone_col = lambda wildcards: name_map[wildcards.cloneShift_method]
    shell: "papermill -p indir {params.indir} -p outdir {params.outdir} -p condition {params.cond} -p donor {params.donor} -p clone_col {params.clone_col} {params.script} {output.note}"


def get_compare_output(wildcards):
    w = wildcards
    if len(condition) > 1:
        return f"{w.outdir}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/{w.cloneShift_method}/knn_kparam_{w.kparam}/output_compare.ipynb"
    else:
        if "inputOnly" in condition:
            return f"{w.outdir}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/{w.cloneShift_method}/knn_kparam_{w.kparam}/output_single_inputOnly_compare.ipynb"
        elif "noInput" in condition:
            return f"{w.outdir}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/{w.cloneShift_method}/knn_kparam_{w.kparam}/output_single_noInput_compare.ipynb"
    return

rule mv_tmp_compare:
    input:
        get_compare_output
    output:
        note = "{outdir}/clonal_shifts/variants_{variants}/donors/donor{d}/{cloneShift_method}/knn_kparam_{kparam}/output.ipynb"
    shell: "cp {input} {output}"


rule compare_input_and_culture_mt:
    input:
        noInput = "{outdir}/clonal_shifts/variants_{variants}/{cloneShift_method}/results/noInput/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
        input = "{outdir}/clonal_shifts/variants_{variants}/{cloneShift_method}/results/inputOnly/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
    output:
        note = "{outdir}/clonal_shifts/variants_{variants}/donors/donor{d}/{cloneShift_method}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output_compare.ipynb"
    params:
        noInput_indir= lambda wildcards, input: dirname(input.noInput),
        input_indir = lambda wildcards, input: dirname(input.input),
        outdir = lambda wildcards, output: dirname(output.note),
        donor = lambda wildcards: wildcards.d,
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/combine_conditions_hypergeometric.ipynb"), #get_script,
        clone_col = lambda wildcards: name_map[wildcards.cloneShift_method]
    shell: "papermill -p noIn_indir {params.noInput_indir} -p input_indir {params.input_indir} -p donor {params.donor} -p clone_col {params.clone_col} -p outdir {params.outdir} {params.script} {output.note}"


rule single_input_and_culture_mt:
    input:
        cond = "{outdir}/clonal_shifts/variants_{variants}/{cloneShift_method}/results/{condition}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
    output:
        note = "{outdir}/clonal_shifts/variants_{variants}/donors/donor{d}/{cloneShift_method}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output_single_{condition}_compare.ipynb"
    params:
        indir= lambda wildcards, input: dirname(input.cond),
        outdir = lambda wildcards, output: dirname(output.note),
        donor = lambda wildcards: wildcards.d,
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/individual_conditions_hypergeometric.ipynb"), #get_script,
        cond = condition[0], # inputOnly or noInput
        clone_col = lambda wildcards: name_map[wildcards.cloneShift_method]
    shell: "papermill -p indir {params.indir} -p outdir {params.outdir} -p condition {params.cond} -p donor {params.donor} -p clone_col {params.clone_col} {params.script} {output.note}"


def get_compare_output_mt(wildcards):
    w = wildcards
    if len(condition) > 1:
        return f"{w.outdir}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/{w.cloneShift_method}/af.{w.af}_othaf.{w.othaf}_cov.{w.cov}_othcov.{w.othcov}_ncells.{w.ncells}_othncells.{w.othncells}_mean.{w.mean}/output_compare.ipynb"
    else:
        if "inputOnly" in condition:
            return f"{w.outdir}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/{w.cloneShift_method}/af.{w.af}_othaf.{w.othaf}_cov.{w.cov}_othcov.{w.othcov}_ncells.{w.ncells}_othncells.{w.othncells}_mean.{w.mean}/output_single_inputOnly_compare.ipynb"
        elif "noInput" in condition:
            return f"{w.outdir}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/{w.cloneShift_method}/af.{w.af}_othaf.{w.othaf}_cov.{w.cov}_othcov.{w.othcov}_ncells.{w.ncells}_othncells.{w.othncells}_mean.{w.mean}/output_single_noInput_compare.ipynb"
    return

rule mv_tmp_compare_mt:
    input:
        get_compare_output_mt
    output:
        note = "{outdir}/clonal_shifts/variants_{variants}/donors/donor{d}/{cloneShift_method}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
    shell: "cp  {input} {output}"



best_p = config["mt_as_clones"]["best_params"]
params_clones = config["clones"]

#print("dendro_d", dendro_d)

rule finalize:
    input:
         expand("{{outdir}}/clonal_shifts/variants_{{variants}}/dendro_bc/results/{condition}/knn/kparam_{kparam}/output.ipynb",
                condition=condition, kparam=params_clones["knn"]["params"]["resolution"]),
         expand("{{outdir}}/clonal_shifts/variants_{{variants}}/clones/results/{condition}/knn/kparam_{kparam}/output.ipynb",
                condition=condition, kparam=params_clones["knn"]["params"]["resolution"]),
         # expand("{{outdir}}/clonal_shifts/variants_{{variants}}/mt_as_clones_dendro/results/{condition}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
         #        condition=condition,
         #        af= best_p["af"], othaf= best_p["oth_af"],
         #        cov= best_p["cov"], othcov= best_p["oth_cov"],
         #        ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
         #        mean= best_p["mean_pos_cov"]),
         # expand("{{outdir}}/clonal_shifts/variants_{{variants}}/mt_as_clones/results/{condition}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
         #        af= best_p["af"], othaf= best_p["oth_af"],
         #        cov= best_p["cov"], othcov= best_p["oth_cov"],
         #        ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
         #        mean= best_p["mean_pos_cov"], condition=condition),
    output:
          "{outdir}/clonal_shifts/variants_{variants}/.complete.txt"
    shell: "touch {output}"




rule finalize_compare_input_and_culture:
    input:
         expand("{{outdir}}/clonal_shifts/variants_{{variants}}/donors/donor{d}/clones/knn_kparam_{kparam}/output.ipynb",
                kparam=params_clones["knn"]["params"]["resolution"], d = np.arange(config["N_DONORS"])),
         expand("{{outdir}}/clonal_shifts/variants_{{variants}}/donors/donor{d}/dendro_bc/knn_kparam_{kparam}/output.ipynb",
                 kparam=params_clones["knn"]["params"]["resolution"], d = np.arange(config["N_DONORS"])), #dt=dendro_d,
         # expand("{{outdir}}/clonal_shifts/variants_{{variants}}/donors/donor{d}/mt_as_clones/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
         #        af= best_p["af"], othaf= best_p["oth_af"],
         #        cov= best_p["cov"], othcov= best_p["oth_cov"],
         #        ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
         #        mean= best_p["mean_pos_cov"], d=np.arange(config["N_DONORS"])),
         # expand("{{outdir}}/clonal_shifts/variants_{{variants}}/donors/donor{d}/mt_as_clones_dendro/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
         #        af= best_p["af"], othaf= best_p["oth_af"],
         #        cov= best_p["cov"], othcov= best_p["oth_cov"],
         #        ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
         #        mean= best_p["mean_pos_cov"], d=np.arange(config["N_DONORS"])),


         # expand("{{outdir}}/clonal_shifts/variants_{{variants}}/dendro_bc/results/{condition}/knn/kparam_{kparam}/dendro_dt_{dt}/output.ipynb",
         #        dt=dendro_d,condition=condition, kparam=params_clones["knn"]["params"]["resolution"]),
         # expand("{{outdir}}/clonal_shifts/variants_{{variants}}/clones/results/{condition}/knn/kparam_{kparam}/output.ipynb",
         #        condition=condition, kparam=params_clones["knn"]["params"]["resolution"]),
         # expand("{{outdir}}/clonal_shifts/variants_{{variants}}/mt_as_clones_dendro/results/{condition}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
         #        condition=condition,
         #        af= best_p["af"], othaf= best_p["oth_af"],
         #        cov= best_p["cov"], othcov= best_p["oth_cov"],
         #        ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
         #        mean= best_p["mean_pos_cov"]),
         # expand("{{outdir}}/clonal_shifts/variants_{{variants}}/mt_as_clones/results/{condition}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
         #        af= best_p["af"], othaf= best_p["oth_af"],
         #        cov= best_p["cov"], othcov= best_p["oth_cov"],
         #        ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
         #        mean= best_p["mean_pos_cov"], condition=condition),
    output:
          "{outdir}/clonal_shifts/variants_{variants}/.finalize_complete.txt"
    shell: "touch {output}"

