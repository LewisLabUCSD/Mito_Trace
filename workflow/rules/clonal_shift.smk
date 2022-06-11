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
# rule all:
#     input:
#         expand("{outdir}/clonal_shifts/{clone_type}/{condition}/sepDon/donor_{d}/results/output.ipynb",
#                outdir = config["outdir"], condition=["inputOnly", "noInput"],
#                clone_type=["clones", dendro_d, "mt_bin", "mt_den"],
#                d=np.arange(config["N_DONORS"])),
#         expand("{outdir}/clonal_shifts/{clone_type}/{condition}/combDon/results/output.ipynb",
#                outdir = config["outdir"], condition=["inputOnly", "noInput"],
#                clone_type=["clones", dendro_d, "mt_bin", "mt_den"])
#
#
# def save_groups_input(df, out_f, clone_col, atac_col, min_clone_size=10):
#     # Get combDonors
#     name_cond_size = df.groupby([clone_col]).size()
#     name_cond_size = name_cond_size[name_cond_size>min_clone_size]
#     clones_filt = name_cond_size.index
#
#     sizes = cells_meta.groupby(clone_col).size().sort_values(ascending=False)
#     sizes = sizes.loc[clones_filt].sort_values(ascending=False)
#     groups = df.groupby([atac_col, clone_col]).size().reset_index().rename({0:"count"}, axis=1)
#     clones = clones_filt#np.unique(groups["name"])
#     atac_cl = np.unique(groups[atac_col])
#     groups.to_csv(out_f,sep="\t")
#     pickle.dump([sizes, clones, atac_cl, clone_col, atac_col], open(str(out_f).replace("groups.tsv", "input_params.p"), "wb"))
#     return
#
#
# rule preprocess_dendro_bc:
#     input:
#         se_meta = "{outdir}/annotation_clones/se_cells_meta.tsv",
#         barcodes_dir = "{outdir}/barcodes/btwnClones/",
#     output:
#         sep = expand("{{outdir}}/clonal_shifts/dendro_{dt}/{condition}/sepDon/donor_{d}/groups.tsv",
#                      d=np.arange(config["N_DONORS"])),
#         combined =  "{outdir}/clonal_shifts/dendro_{dt}/{condition}/combDon/groups.tsv"
#     params:
#         N_DONORS = config["N_DONORS"],
#         min_clone_size = 10
#     run:
#         clone_col = "den_clust"
#         barcodes_in = {}
#         for d in np.arange(params.N_DONORS):
#             barcodes_in[d] = pd.read_csv(join(input.barcodes_dir,f"donor{d}.clones_dendro.csv"), index_col=0)
#
#             barcodes_in[d][clone_col] = str(d) + "_" + barcodes_in[d][clone_col]
#
#         cells_meta = pd.read_csv(input.se_cells_meta_f, sep="\t")
#         cells_meta = cells_meta.loc[~(cells_meta["name"]=="None")]
#
#         if not "cluster_labels" in cells_meta.columns.values:
#             cells_meta["cluster_labels"] = cells_meta["seurat_clusters"]
#
#         cells_meta[clone_col] = cells_meta.apply(lambda x: barcodes_in[int(x["donor"])].loc[x["name"], clone_col] , axis=1)
#
#         if wildcards.condition == "inputOnly":
#             cells_meta = cells_meta.loc[cells_meta["condition"]==params.input_cond]
#         elif wildcards.condition == "noInput":
#             cells_meta = cells_meta.loc[~(cells_meta["condition"]==params.input_cond)]
#         save_groups_input(cells_meta, output.combined, clone_col, "cluster_labels", min_clone_size=params.min_clone_size)
#         # Save sepDonors
#         for d, val in cells_meta.groupby("donor"):
#             print(d)
#             save_groups_input(val, output.sep[int(d)], clone_col, "cluster_labels", min_clone_size=params.min_clone_size)
#
#
#
# rule preprocess_clones:
#     input: "{outdir}/annotation_clones/se_cells_meta.tsv"
#     output:
#         sep = expand("{{outdir}}/clonal_shifts/clones/{condition}/sepDon/donor_{d}/groups.tsv",
#                      d=np.arange(config["N_DONORS"])),
#         combined =  "{outdir}/clonal_shifts/clones/{condition}/combDon/groups.tsv"
#     params:
#         input_cond = "Input",
#         min_clone_size = 10
#     run:
#         cells_meta = pd.read_csv(input.se_cells_meta_f, sep="\t")
#         cells_meta = cells_meta.loc[~(cells_meta["name"]=="None")]
#
#         if not "cluster_labels" in cells_meta.columns.values:
#             cells_meta["cluster_labels"] = cells_meta["seurat_clusters"]
#
#         sizes = cells_meta.groupby("name").size().sort_values(ascending=False)
#         if wildcards.condition == "inputOnly":
#             cells_meta = cells_meta.loc[cells_meta["condition"]==params.input_cond]
#         elif wildcards.condition == "noInput":
#             cells_meta = cells_meta.loc[~(cells_meta["condition"]==params.input_cond)]
#
#         # Save combDonors
#         save_groups_input(cells_meta, output.combined, "name", "cluster_labels", min_clone_size=params.min_clone_size)
#
#         # Save sepDonors
#         for d, val in cells_meta.groupby("donor"):
#             print(d)
#             save_groups_input(val, output.sep[int(d)], "name", "cluster_labels", min_clone_size=params.min_clone_size)
#
#
# rule preprocess_mt_bin:
#     input: "{outdir}/mt_clones_thresh/best_params/best_params_save.ipynb"
#     output: "{outdir}/clonal_shifts/mt_bin/groups.tsv"
#     params:
#         indir = lambda wildcards, input: dirname(input[0]),
#     run:
#         df = pd.read_csv(params.indir)
#
# #mt_clones_thresh/dendro
# rule preprocess_mt_dendro:
#     input: "{outdir}/mt_clones_thresh/dendro/don_{d}_mt_dendro_clust.csv"
#     output:
#         sep = expand("{{outdir}}/clonal_shifts/mt_den/{condition}/sepDon/donor_{d}/groups.tsv",
#                      d=np.arange(config["N_DONORS"])),
#         combined =  "{outdir}/clonal_shifts/mt_den/{condition}/combDon/groups.tsv"
#         #"{outdir}/clonal_shifts/mt_den/groups.tsv"
#     shell: ""
#
#
# rule clonal_shift:
#     input: "{outdir}/groups.tsv"
#     output: "{outdir}/results/output.ipynb"
#     params:
#         outdir = lambda wildcards, output: dirname(output[0]),
#         note = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/run_hypergeo.ipynb"),
#     shell: "papermill -p groups {input} -p outdir {params.outdir} {params.note} {output.note} "
#
#
# rule run_sep:
#     input:
#         "{outdir}/clonal_shifts/{clone_type}/{condition}/sepDon/donor_{d}/groups.tsv"
#     output:
#         "{outdir}/clonal_shifts/{clone_type}/{condition}/sepDon/donor_{d}/results/output.ipynb"
#     params:
#         outdir = lambda wildcards, output: dirname(output[0]),
#         note = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/run_hypergeo.ipynb")
#     shell: "papermill -p groups {input} -p outdir {params.outdir} {params.note} {output.note} "
#
#
# rule run_comb:
#     input:
#         "{outdir}/clonal_shifts/{clone_type}/{condition}/combDon/groups.tsv"
#     output:
#         "{outdir}/clonal_shifts/{clone_type}/{condition}/combDon/results/output.ipynb"
#     params:
#         outdir = lambda wildcards, output: dirname(output[0]),
#         note = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/run_hypergeo.ipynb")
#     shell: "papermill -p groups {input} -p outdir {params.outdir} {params.note} {output.note} "
#
#
#
# rule finalize:
#     input:
#         expand("{outdir}/clonal_shifts/{clone_type}/{condition}/sepDon/donor_{d}/results/output.ipynb",
#                outdir = config["outdir"], condition=["inputOnly", "noInput"],
#                clone_type= ["clones"],#["clones", dendro_d, "mt_bin", "mt_den"],
#                d=np.arange(config["N_DONORS"])),
#         expand("{outdir}/clonal_shifts/{clone_type}/{condition}/combDon/results/output.ipynb",
#                outdir = config["outdir"], condition=["inputOnly", "noInput"],
#                clone_type= ["clones"])#["clones", dendro_d, "mt_bin", "mt_den"],
#

# def get_script(wildcards):
#     clone_type = wildcards.clone_type
#     if clone_type == "clones":
#         return join(ROOT_DIR, "workflow/notebooks/clonal_shifts/hypergeometric_clones.ipynb")
#     elif clone_type == "mt_den":
#         return join(ROOT_DIR, "workflow/notebooks/clonal_shifts/hypergeometric_mt_as_clones_dendro.ipynb")
#     elif clone_type == "mt_bin":
#         return join(ROOT_DIR, "workflow/notebooks/clonal_shifts/hypergeometric_mt_as_clones.ipynb")
#     elif clone_type == "dendro_bc":
#         return join(ROOT_DIR, "workflow/notebooks/clonal_shifts/dendro.ipynb")
#     raise ValueError("clone_type variable")
#     return


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
        "{outdir}/clonal_shifts/variants_{variants}/clones/results/{condition}/knn/kparam_{kparam}/output.ipynb"
    params:
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/hypergeometric_clones.ipynb"), #get_script,
        outdir = lambda wildcards, output: dirname(output[0]),
        N_DONORS = config["N_DONORS"],
    shell: "papermill -p se_cells_meta_f {input.se_meta} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} -p condition {wildcards.condition} {params.script} {output}"

rule compare_input_and_culture_clones:
    input:
        noInput = "{outdir}/clonal_shifts/variants_{variants}/clones/results/noInput/knn/kparam_{kparam}/output.ipynb",
        input = "{outdir}/clonal_shifts/variants_{variants}/clones/results/inputOnly/knn/kparam_{kparam}/output.ipynb"
    output:
        note = "{outdir}/clonal_shifts/variants_{variants}/donors/donor{d}/clones/knn_kparam_{kparam}/output.ipynb"
    params:
        noInput_indir= lambda wildcards, input: dirname(input.noInput),
        input_indir = lambda wildcards, input: dirname(input.input),
        outdir = lambda wildcards, output: dirname(output.note),
        donor = lambda wildcards: wildcards.d,
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/combine_conditions_hypergeometric.ipynb"), #get_script,
    shell: "papermill -p noIn_indir {params.noInput_indir} -p input_indir {params.input_indir} -p donor {params.donor} -p clone_col name -p outdir {params.outdir} {params.script} {output.note}"


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


rule compare_input_and_culture_dendro_bc:
    input:
        noInput = "{outdir}/clonal_shifts/variants_{variants}/dendro_bc/results/noInput/knn/kparam_{kparam}/output.ipynb",
        input = "{outdir}/clonal_shifts/variants_{variants}/dendro_bc/results/inputOnly/knn/kparam_{kparam}/output.ipynb"
    output:
        note = "{outdir}/clonal_shifts/variants_{variants}/donors/donor{d}/dendro_bc/knn_kparam_{kparam}/output.ipynb"
    params:
        noInput_indir= lambda wildcards, input: dirname(input.noInput),
        input_indir = lambda wildcards, input: dirname(input.input),
        outdir = lambda wildcards, output: dirname(output.note),
        donor = lambda wildcards: wildcards.d,
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/combine_conditions_hypergeometric.ipynb"), #get_script,
    shell: "papermill -p noIn_indir {params.noInput_indir} -p input_indir {params.input_indir} -p donor {params.donor} -p clone_col den_clust -p outdir {params.outdir} {params.script} {output.note}"



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


rule compare_input_and_culture_mt_as_clones:
    input:
        noInput = "{outdir}/clonal_shifts/variants_{variants}/mt_as_clones/results/noInput/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
        input = "{outdir}/clonal_shifts/variants_{variants}/mt_as_clones/results/inputOnly/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
    output:
        note = "{outdir}/clonal_shifts/variants_{variants}/donors/donor{d}/mt_as_clones/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
    params:
        noInput_indir= lambda wildcards, input: dirname(input.noInput),
        input_indir = lambda wildcards, input: dirname(input.input),
        outdir = lambda wildcards, output: dirname(output.note),
        donor = lambda wildcards: wildcards.d,
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/combine_conditions_hypergeometric.ipynb"), #get_script,
    shell: "papermill -p noIn_indir {params.noInput_indir} -p input_indir {params.input_indir} -p donor {params.donor} -p clone_col Variants -p outdir {params.outdir} {params.script} {output.note}"


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


rule compare_input_and_culture_mt_as_clones_dendro:
    input:
        noInput = "{outdir}/clonal_shifts/variants_{variants}/mt_as_clones_dendro/results/noInput/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
        input = "{outdir}/clonal_shifts/variants_{variants}/mt_as_clones_dendro/results/inputOnly/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
    output:
        note = "{outdir}/clonal_shifts/variants_{variants}/donors/donor{d}/mt_as_clones_dendro/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb"
    params:
        noInput_indir= lambda wildcards, input: dirname(input.noInput),
        input_indir = lambda wildcards, input: dirname(input.input),
        outdir = lambda wildcards, output: dirname(output.note),
        donor = lambda wildcards: wildcards.d,
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/combine_conditions_hypergeometric.ipynb"), #get_script,
    shell: "papermill -p noIn_indir {params.noInput_indir} -p input_indir {params.input_indir} -p donor {params.donor} -p clone_col den_clust -p outdir {params.outdir} {params.script} {output.note}"



best_p = config["mt_as_clones"]["best_params"]
params_clones = config["clones"]

print("dendro_d", dendro_d)

rule finalize:
    input:
         expand("{{outdir}}/clonal_shifts/variants_{{variants}}/dendro_bc/results/{condition}/knn/kparam_{kparam}/output.ipynb",
                condition=condition, kparam=params_clones["knn"]["params"]["resolution"]),
         expand("{{outdir}}/clonal_shifts/variants_{{variants}}/clones/results/{condition}/knn/kparam_{kparam}/output.ipynb",
                condition=condition, kparam=params_clones["knn"]["params"]["resolution"]),
         expand("{{outdir}}/clonal_shifts/variants_{{variants}}/mt_as_clones_dendro/results/{condition}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
                condition=condition,
                af= best_p["af"], othaf= best_p["oth_af"],
                cov= best_p["cov"], othcov= best_p["oth_cov"],
                ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
                mean= best_p["mean_pos_cov"]),
         expand("{{outdir}}/clonal_shifts/variants_{{variants}}/mt_as_clones/results/{condition}/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
                af= best_p["af"], othaf= best_p["oth_af"],
                cov= best_p["cov"], othcov= best_p["oth_cov"],
                ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
                mean= best_p["mean_pos_cov"], condition=condition),
    output:
          "{outdir}/clonal_shifts/variants_{variants}/.complete.txt"
    shell: "touch {output}"



rule finalize_compare_input_and_culture:
    input:
         expand("{{outdir}}/clonal_shifts/variants_{{variants}}/donors/donor{d}/clones/knn_kparam_{kparam}/output.ipynb",
                kparam=params_clones["knn"]["params"]["resolution"], d = np.arange(config["N_DONORS"])),
         expand("{{outdir}}/clonal_shifts/variants_{{variants}}/donors/donor{d}/dendro_bc/knn_kparam_{kparam}/output.ipynb",
                 kparam=params_clones["knn"]["params"]["resolution"], d = np.arange(config["N_DONORS"])), #dt=dendro_d,
         expand("{{outdir}}/clonal_shifts/variants_{{variants}}/donors/donor{d}/mt_as_clones/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
                af= best_p["af"], othaf= best_p["oth_af"],
                cov= best_p["cov"], othcov= best_p["oth_cov"],
                ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
                mean= best_p["mean_pos_cov"], d=np.arange(config["N_DONORS"])),
         expand("{{outdir}}/clonal_shifts/variants_{{variants}}/donors/donor{d}/mt_as_clones_dendro/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/output.ipynb",
                af= best_p["af"], othaf= best_p["oth_af"],
                cov= best_p["cov"], othcov= best_p["oth_cov"],
                ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
                mean= best_p["mean_pos_cov"], d=np.arange(config["N_DONORS"])),
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

