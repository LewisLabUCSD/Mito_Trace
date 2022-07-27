from src.config import ROOT_DIR
from os.path import join, dirname
import numpy as np
import pandas as pd

# wildcard_constraints:
#     outdir='.+enriched_barcodes$'

knn_clone_shift= ["clones"] #knn_clone_shift= ["clones", "dendro_bc"]

params_clones = config["clones"]
sample_ids = config["samples"].index
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

gff = config["genome_path"][config["genome"]]["gff"]
dendro_d = params_clones['params']['dendro_thresh']

objectives_l = params_clones["params"]["clone_barcodes_optimization"]["objectives"]
weights_cfg = params_clones["params"]["clone_barcodes_optimization"]["weights_cfg"]

#
# def get_anno_cells_meta(wildcards):
#     w = wildcards
#     return f"{w.outdir}/gff_{gff}/annotation_clones/se_cells_meta_labels.tsv"
#
import os
def get_anno_cells(wildcards):
    w = wildcards
    return  f"{w.outdir}/clones/variants_{w.variants}/knn/kparam_{w.kparam}/gff_{gff}/annotation_clones/se_cells_meta_labels.tsv"


rule run_variants_params:
    """ Run the parameter optimization rule.
    
        Get the best result and output the cells meta, af and dp.
        For now also requires the labelled cells_meta file 
        for filtering later.
    """
    input:
        clone_cells = "{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
        anno_cells_meta_f = get_anno_cells #expand("{{outdir}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/gff_{gff}/annotation_clones/se_cells_meta_labels.tsv",
                                   #gff=gff) # /cells_meta_labels.tsv"
    output:
        #params_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/params.csv",
        #all_params_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/raw_params.csv",
        params_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/params.csv",
        objs_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objectives.csv"
        # best_params_f ="{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/best_params.csv",
        # best_params_clone_vars_f ="{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/best_params_filt_clone_vars.csv",
        # cells_meta_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
        # af_f ="{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/af.tsv",
        # dp_f ="{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/dp.tsv",
    params:
        outdir = lambda wildcards, output: dirname(output.params_f),
        indir = lambda wildcards, input: dirname(input.clone_cells),
        donor = lambda wildcards: wildcards.d,
        # Objective weights. order of the columns
        weights = [1,1,1,1,-1, 1, 1], #weights=lambda wildcards: config["enriched_vars"]["params"]["objective_weights"][wildcards.objective_id]
        objectives_l = config[""] #["variants_with_clone_norm_by_1_over_nclones_with_variant", "max_clone_ncells_over_ncells","pct_thresh","other_pct_thresh", "n_vars", "obj_nclones_more_than_one_unique"],
        ncpus=12, #config["ncpus"]
        topn=16,
        to_test=False
    threads: 17
    priority: 8
    script: join(ROOT_DIR, "workflow/notebooks/clone_vars/optimization_run.py.py")
    # log:
    #     notebook = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/output.ipynb",
    # notebook:
    #     join(ROOT_DIR, "workflow/notebooks/clone_vars/optimization_run.py.ipynb")


# rule get_weights:
#     input:
#         clone_cells = "{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
#         all_params_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/raw_params.csv",
#         objs_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objectives.csv"
#     output:
#         params_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/params.csv",
#         loss_multi_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/loss_multi.pdf",
#         top_params_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/top_param_results.pdf",
#     params:
#         indir = lambda wildcards, input: dirname(input.clone_cells),
#         outdir = lambda wildcards, output: dirname(output.params_f),
#         donor = lambda wildcards: wildcards.d,
#         weights = [1,1,1,1,-1, 1, 1], #weights=lambda wildcards: config["enriched_vars"]["params"]["objective_weights"][wildcards.objective_id]
#         objectives_l = ["variants_with_clone_norm_by_1_over_nclones_with_variant", "max_clone_ncells_over_ncells",
#                         "pct_thresh","other_pct_thresh", "n_vars", "obj_nclones_more_than_one_unique"],
#         topn=16,
#         to_test=False
#     priority: 10
#     log: "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/weights.ipynb",
#     notebook:
#         join(ROOT_DIR, "workflow/notebooks/clone_vars/weights_of_objectives.py.ipynb")
#

rule optim_results:
    input:
        params_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/params.csv",
        anno_cells_meta_f = get_anno_cells,
        clone_cells = "{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
    output:
        #note = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/out.ipynb",
        best_params_f ="{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/best_params.csv",
        best_params_clone_vars_f ="{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/best_params_filt_clone_vars.csv",
        cells_meta_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
        af_f ="{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/af.tsv",
        dp_f ="{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/dp.tsv",
        top_clustered = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/top_param_group_results.pdf",
        top_table = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/top_param_groups_clone_vars.pdf",
        top_clustered_table = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/top10perc_param_group_clusters.pdf",
    priority: 10
    params:
        indir = lambda wildcards, input: dirname(input.params_f),
        outdir = lambda wildcards, output: dirname(output.top_clustered),
        af_indirs = lambda wildcards, input: dirname(input.clone_cells),
        #note = join(ROOT_DIR, "workflow/notebooks/clone_vars/optimization_results.ipynb"),
        weights = [1,1,1,1,-1, 1, 1],
        objectives_l = ["variants_with_clone_norm_by_1_over_nclones_with_variant", "max_clone_ncells_over_ncells",
                        "pct_thresh","other_pct_thresh", "n_vars", "obj_nclones_more_than_one_unique"],
    log:
        notebook = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/out.ipynb"
    notebook:
        join(ROOT_DIR, "workflow/notebooks/clone_vars/multi_cluster_optimization_results.py.ipynb")
        #join(ROOT_DIR, "workflow/notebooks/clone_vars/optimization_results.py.ipynb")


rule merge_donor_optim:
    input:
        cells_meta_f = expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{d}/cells_meta.tsv", d=np.arange(config["N_DONORS"])),
        af_f =expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{d}/af.tsv", d=np.arange(config["N_DONORS"])),
        dp_f =expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{d}/dp.tsv", d=np.arange(config["N_DONORS"])),
    output:
        cells_meta_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
        af_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/af.tsv",
        dp_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/dp.tsv",
    run:
        pd.concat([pd.read_csv(x, sep="\t", index_col=0) for x in input.cells_meta_f]).to_csv(output.cells_meta_f, sep="\t")
        pd.concat([pd.read_csv(x, sep="\t", index_col=0) for x in input.af_f], axis=0).to_csv(output.af_f, sep="\t")
        pd.concat([pd.read_csv(x, sep="\t", index_col=0) for x in input.dp_f], axis=0).to_csv(output.dp_f, sep="\t")


rule barcodes_btwnClones_dendro:
    input:
        cells_meta = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
        af_note = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/af.tsv",
    output:
        note = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.ipynb",
        mean = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.mean.csv",
        res = report(multiext("{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.",
                                  "clones_dendro.csv", "dendrogram_pvals.txt",
                                  "dendro.NoCondition.max2.AF.png"),
                         category="lineage")
    params:
        note = join(ROOT_DIR, "workflow", "notebooks", "clone_af_dendrograms", "MT_btwnClones_Barcode_dendro_dynamicClust.ipynb"),
        indir = lambda wildcards, input: dirname(input.cells_meta),
        outdir = lambda wildcards, output: dirname(output.note),
        #dendro_thresh = 0.6
    shell: "papermill -p INDIR {params.indir} -p OUTDIR {params.outdir} -p DONOR {wildcards.d} -p dendroThresh {wildcards.dendro_thresh} {params.note} {output.note}"

# rule clonalshift_clones:
#     input:
#         se_meta = expand("{{outdir}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/gff_{gff}/annotation_clones/se_cells_meta_labels.tsv",
#                          gff=config["gff"])
#     output:
#         note = "{outdir}/clonal_shifts/variants_{variants}/clones/results/{condition}/knn/kparam_{kparam}/output.ipynb"
#     params:
#         script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/hypergeometric_clones.ipynb"), #get_script,
#         outdir = lambda wildcards, output: dirname(output.note),
#         N_DONORS = config["N_DONORS"],
#     shell: "papermill -p se_cells_meta_f {input.se_meta} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} -p condition {wildcards.condition} {params.script} {output}"
#

name_map = {"clones":"name", "dendro_bc": "den_clust", "mt_as_clones": "Variants", "mt_as_clones_dendro": "den_clust"}

module cloneShiftMod:
    snakefile: "../clonal_shift.smk"
    config: config

#use rule * from cloneShiftMod as cloneshifttwo_*

module indClonesMod:
    snakefile: "../individual_clones.smk"
    config: config

#use rule * from indClonesMod as indclonestwo_*


rule counts_clones:
    input:
        #se_meta = "{outdir}/annotation_clones/se_cells_meta.tsv",
        se_meta = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
    output:
        multiext("{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/annotation_clones/clone_counts/",
                 "clone_counts.barplot_conditions.png",
                 ),
        note = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/annotation_clones/clone_counts/counts_clones.ipynb"
    params:
        outdir = lambda wildcards, output: dirname(output.note),
        note = join(ROOT_DIR, "workflow/notebooks/lineage_clones/python_clone_cell_counts.ipynb"),
        sample_names = ",".join(config["samples"].index),
        min_cell = config["annotation_clones"]["params"]["clone_sizes_min_cell"]
    shell:
        "papermill -p se_cells_meta_f {input.se_meta} -p outdir {params.outdir} -p sample_names {params.sample_names} -p min_cell {params.min_cell} {params.note} {output.note}"


## Clonal shift rules:
## a) rule clonalshift_clones  and clonalshift_dendro_bc
## b) rule finalize
#######################################################################
use rule clonalshift_clones from cloneShiftMod as cloneshifttwo_clonalshift_clones with:
    input:
        se_meta = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv", #"{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
        # se_meta = expand("{{outdir}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/annotation_clones/se_cells_meta_labels.tsv",
        #              gff=config["gff"])
    output:
        note = "{outdir}/enriched_barcodes/clonal_shifts/variants_{variants}/clones/results/{condition}/knn/kparam_{kparam}/output.ipynb"

# use rule clonalshift_clonalshift_dendro_bc from cloneShiftMod as cloneshifttwo_clonalshift_dendro_bc with:
#     input:
#         se_meta = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
#         barcodes_dir = expand("{{outdir}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/barcodes/btwnClones_dendro_dt_{dt}/donor{d}.clones_dendro.csv",
#                               dt=dendro_d, d=np.arange(config["N_DONORS"]))
#     output:
#         note = "{outdir}/enriched_barcodes/clonal_shifts/variants_{variants}/dendro_bc/results/{condition}/knn/kparam_{kparam}/output.ipynb"
#

use rule compare_input_and_culture_cl from cloneShiftMod as cloneshifttwo_compare_input_and_culture_cl with:
    input:
        noInput = "{outdir}/enriched_barcodes/clonal_shifts/variants_{variants}/clones/results/noInput/knn/kparam_{kparam}/output.ipynb",
        input = "{outdir}/enriched_barcodes/clonal_shifts/variants_{variants}/clones/results/inputOnly/knn/kparam_{kparam}/output.ipynb"
    output:
        note = "{outdir}/enriched_barcodes/clonal_shifts/variants_{variants}/donors/donor{d}/clones/knn_kparam_{kparam}/output_compare.ipynb"
    params:
        clone_col = name_map["clones"], #lambda wildcards: name_map[wildcards.cloneShift_method]
        noInput_indir= lambda wildcards, input: dirname(input.noInput),
        input_indir = lambda wildcards, input: dirname(input.input),
        outdir = lambda wildcards, output: dirname(output.note),
        donor = lambda wildcards: wildcards.d,
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/combine_conditions_hypergeometric.ipynb"), #get_script,



use rule single_input_and_culture_cl from cloneShiftMod as cloneshifttwo_single_input_and_culture_cl with:
    input:
        cond = "{outdir}/enriched_barcodes/clonal_shifts/variants_{variants}/clones/results/{condition}/knn/kparam_{kparam}/output.ipynb",
    output:
        note = "{outdir}/enriched_barcodes/clonal_shifts/variants_{variants}/donors/donor{d}/clones/knn_kparam_{kparam}/output_single_{condition}_compare.ipynb"
    params:
        indir= lambda wildcards, input: dirname(input.cond),
        outdir = lambda wildcards, output: dirname(output.note),
        donor = lambda wildcards: wildcards.d,
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/individual_conditions_hypergeometric.ipynb"), #get_script,
        cond = condition[0],
        clone_col = name_map["clones"], #lambda wildcards: name_map[wildcards.cloneShift_method]






def get_compare_output(wildcards):
    w = wildcards
    if len(condition) > 1:
        return f"{w.outdir}/enriched_barcodes/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/clones/knn_kparam_{w.kparam}/output_compare.ipynb"
    else:
        if "inputOnly" in condition:
            return f"{w.outdir}/enriched_barcodes/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/clones/knn_kparam_{w.kparam}/output_single_inputOnly_compare.ipynb"
        elif "noInput" in condition:
            return f"{w.outdir}/enriched_barcodes/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/clones/knn_kparam_{w.kparam}/output_single_noInput_compare.ipynb"
    return

use rule mv_tmp_compare from cloneShiftMod as cloneshifttwo_mv_tmp_compare with:
    input:
        get_compare_output
    output:
        note = "{outdir}/enriched_barcodes/clonal_shifts/variants_{variants}/donors/donor{d}/clones/knn_kparam_{kparam}/output.ipynb"

# use rule finalize from cloneShiftMod as cloneshifttwo_finalize with:
#     input:
#          expand("{{outdir}}/enriched_barcodes/clonal_shifts/enriched_barcodes/variants_{{variants}}/dendro_bc/results/{condition}/knn/kparam_{kparam}/output.ipynb",
#                 condition=condition, kparam=params_clones["knn"]["params"]["resolution"]),
#          expand("{{outdir}}/enriched_barcodes/clonal_shifts/variants_{{variants}}/clones/results/{condition}/knn/kparam_{kparam}/output.ipynb",
#                 condition=condition, kparam=params_clones["knn"]["params"]["resolution"]),
#     output:
#           "{outdir}/enriched_barcodes/clonal_shifts/variants_{variants}/.complete.txt"
#     shell: "touch {output}"
#
use rule finalize_compare_input_and_culture from cloneShiftMod as cloneshifttwo_finalize_compare_input_and_culture with:
    input:
         expand("{{outdir}}/enriched_barcodes/clonal_shifts/enriched_barcodes/variants_{{variants}}/dendro_bc/results/{condition}/knn/kparam_{kparam}/output.ipynb",
                condition=condition, kparam=params_clones["knn"]["params"]["resolution"]),
         expand("{{outdir}}/enriched_barcodes/clonal_shifts/variants_{{variants}}/clones/results/{condition}/knn/kparam_{kparam}/output.ipynb",
                condition=condition, kparam=params_clones["knn"]["params"]["resolution"]),
    output:
          "{outdir}/enriched_barcodes/clonal_shifts/variants_{variants}/.finalize_complete.txt"
    #shell: "touch {output}"



## Indclones
## a) rule get_clone_cells and get_clone_dendro_bc
## b) top_clone_complete
#######################################################################
use rule get_clone_cells from indClonesMod as indclonestwo_get_clone_cells with:
    input:
        #se_meta = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/se_cells_meta_labels.tsv"
        se_cells_meta_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv", #"{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
        # se_meta = expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/annotation_clones/cells_meta.tsv",
        #              gff=config["gff"])
    output:
        cells_meta = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv"



# use rule get_clone_dendro_cells from indClonesMod as indclonestwo_get_clone_dendro_cells with:
#     input:
#         barcodes_dir = expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/barcodes/btwnClones_dendro_dt_{dt}/donor{{d}}.clones_dendro.csv", dt=dendro_d),
#         se_cells_meta_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/annotation_clones/cells_meta.tsv"
#     output:
# #         cells = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_dendro_bc/cells_meta.tsv"
# def get_hypergeo(wildcards):
#     w = wildcards
#     clshift = w.clone_shift_method
#     indir = f"{w.outdir}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/{w.clone_shift_method}"
#     #output.ipynb
#     return f"{indir}/knn_kparam_{w.kparam}/output.ipynb"

use rule get_rank_cl_sizes from indClonesMod as indclonestwo_get_rank_cl_sizes with:
    input:
        hyper = "{outdir}/enriched_barcodes/clonal_shifts/variants_{variants}/donors/donor{d}/clones/knn_kparam_{kparam}/output.ipynb", #get_hypergeo,
        cells="{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv"
    output:
        note = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/clones_ranked/output_rank_ncells.ipynb",
        clone_order = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/clones_ranked/cloneID_rank_ncells.txt"



use rule top_clone_mt_variants_cl from indClonesMod as indclonestwo_top_clone_mt_variants_cl with:
    input:
        #se =  expand("{{outdir}}/anno_multiplex/gff_{gff}/SE.rds", gff=config["gff"]),
        cells = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv",
        clone_order_f = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/clones_ranked/cloneID_rank_ncells.txt",
        mt_f =  "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/af.tsv",
    output:
        note = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/{clone_order}_clone_mt_variants.ipynb",


# use rule top_clone_mt_variants_mt from indClonesMod as indclonestwo_clone_mt_variants_mt with:
#     input:
#         se =  expand("{{outdir}}/anno_multiplex/gff_{gff}/SE.rds", gff=config["gff"]),
#         cells = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv",
#         mt_f =  "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/af.tsv",
#     output:
#         note = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_clones/top/{clone_order}_clone_mt_variants.ipynb",
#

use rule top_clone_umap_overlay from indClonesMod as indclonestwo_top_clone_umap_overlay with:
    input:
        se =  expand("{{outdir}}/anno_multiplex/gff_{gff}/SE.rds", gff=config["gff"]),
        cells = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv",
        clone_order_f = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/clones_ranked/cloneID_rank_ncells.txt",
    output:
        note = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/top_umap_overlay.ipynb"
    params:
        script = join(ROOT_DIR, "workflow/notebooks/individual_clones/top_individual_clones_umap_overlap.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),

# use rule ind_clone_umap_overlay from indClonesMod as indclonestwo_ind_clone_umap_overlay with:
#     input:
#         se =  expand("{{outdir}}/anno_multiplex/gff_{gff}/SE.rds", gff=config["gff"])
#         #cells = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv",
#         #clone_order_f = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/clones_ranked/cloneID_rank_ncells.txt",
#

def get_top_clone_hypergeoSig_script(wildcards):
    return join(ROOT_DIR, "workflow/notebooks/individual_clones/oneCondition_top_individual_clone_lineage_Sig_hypergeo.ipynb")

def get_param_top_clone_hypergeoSig_condition():
    if len(condition) == 1:
        return condition[0]
    else:
        return "None"


use rule top_clone_hypergeoSig_cl from indClonesMod as indclonestwo_top_clone_hypergeoSig_cl with:
    input:
        indir = "{outdir}/enriched_barcodes/clonal_shifts/variants_{variants}/donors/donor{d}/clones/knn_kparam_{kparam}/output.ipynb", #get_hypergeo,
        cloneIDs_f = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/clones_ranked/cloneID_rank_ncells.txt",
        #cloneID = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/cloneIDs/{cloneID}.txt"
    output:
        note = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/top/top_hypergeo_sig.ipynb",


use rule top_clone_lineage_count from indClonesMod as indclonestwo_top_clone_lineage_count with:
    input:
        cells = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv",
        clone_order_f = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/clones_ranked/cloneID_rank_ncells.txt",
    output:
        note = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/top_clone_lineage_count.ipynb",


use rule top_merge from indClonesMod as indclonestwo_top_merge with:
    input:
        ins = multiext("{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/","top_umap_overlay.ipynb", "top_hypergeo_sig.ipynb", "top_clone_lineage_count.ipynb", "top_clone_mt_variants.ipynb"),
    output:
        out_f = multiext("{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/","clone_shift_combine.pdf"),
        note = "{outdir}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/merge_panels.ipynb"


use rule top_clone_complete from indClonesMod as indclonestwo_top_clone_complete with:
    input:
        expand("{{outdir}}/enriched_barcodes/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/{f}",
              f=["top_umap_overlay.ipynb", "top_hypergeo_sig.ipynb", "top_clone_lineage_count.ipynb", "all_clone_mt_variants.ipynb", "clone_shift_combine.pdf"],
              d=np.arange(config["N_DONORS"]), variants=[x for x in params_clones["variants"] if x != "simple"],
              kparam=params_clones["knn"]["params"]["resolution"])
    output:
        "{outdir}/enriched_barcodes/single_clones/.top_clones.txt"




rule finalize:
    input:
        "{outdir}/enriched_barcodes/single_clones/.top_clones.txt",
        expand("{{outdir}}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/annotation_clones/clone_counts/counts_clones.ipynb",
               variants=[x for x in params_clones["variants"] if x != "simple"],
               kparam=params_clones["knn"]["params"]["resolution"])
    output: final = "{outdir}/enriched_barcodes/.finalize"
    shell: "touch {output}"


# rule run_get_best_variants_param:
#     input:
#         indir = "{indir}/clones/variants_{variants}/knn/kparam_30/"
#     output:
#         outdir = "{indir}/clones/variants_{variants}/knn/kparam_30/enriched_barcodes/donor{d}/params_results.tsv"
#     params:
#         donor = lambda wildcards: wildcards.d,
#         # Objective weights. order of the columns
#         weights = [1,0,0,1,-1, 1, 1], #weights=lambda wildcards: config["enriched_vars"]["params"]["objective_weights"][wildcards.objective_id]
#         objectives_l = ["variants_with_clone_norm_by_1_over_nclones_with_variant", "max_clone_ncells_over_ncells",
#                         "pct_thresh","other_pct_thresh", "n_vars", "obj_nclones_more_than_one_unique"],
#         ncpus=8, #config["ncpus"]
#         topn=16
#     notebook:
#         join(ROOT_DIR, "/workflow/notebooks/clone_vars/table_distinguishing_vars_multi_a_and_pctthresh.ipynb")
