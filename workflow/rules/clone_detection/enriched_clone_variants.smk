from src.config import ROOT_DIR
from os.path import join, dirname

rule run_variants_params:
    input:
        indir = "{indir}/clones/variants_init/knn/kparam_30/"
    output:
        params_f= "{indir}/clones/variants_init/knn/kparam_30/enriched_variants/donor{d}/params_results.tsv",
        best_params_f ="{indir}/clones/variants_init/knn/kparam_30/enriched_variants/donor{d}/best_params.tsv",
        cells_meta_f ="{indir}/clones/variants_init/knn/kparam_30/enriched_variants/donor{d}/cells_meta.tsv",
        af_f ="{indir}/clones/variants_init/knn/kparam_30/enriched_variants/donor{d}/af.tsv",
        dp_f ="{indir}/clones/variants_init/knn/kparam_30/enriched_variants/donor{d}/dp.tsv",
    params:
        outdir = lambda wildcards, output: dirname(output.params_f),
        donor = lambda wildcards: wildcards.d,
        # Objective weights. order of the columns
        weights = [1,0,0,1,-1, 1, 1], #weights=lambda wildcards: config["enriched_vars"]["params"]["objective_weights"][wildcards.objective_id]
        objectives_l = ["variants_with_clone_norm_by_1_over_nclones_with_variant",
                        "max_clone_ncells_over_nclones", "max_clone_ncells_over_ncells",
                        "pct_thresh","other_pct_thresh", "n_vars", "obj_nclones_more_than_one_unique"],
        ncpus=8, #config["ncpus"]
        topn=16
    notebook:
        join(ROOT_DIR, "/workflow/notebooks/clone_vars/table_distinguishing_vars_multi_a_and_pctthresh.ipynb")



module cloneShiftMod:
    snakefile: "./rules/clonal_shift.smk"
    config: params

use rule * from cloneShiftMod as cloneshift_*

module indClonesMod:
    snakefile: "./rules/individual_clones.smk"
    config: params

# use rule nuclear_and_mtclone_counts from mtcloneMod as mtclone_nuclear_and_mtclone_counts with:
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_{variant}/anno_mt_clones_thresh/anno_mt_clones.ipynb",
#
#
use rule * from indClonesMod as indclones_*

# run individual clones and use these cells meta.
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


# rule run_get_best_variants_param:
#     input:
#         indir = "{indir}/clones/variants_init/knn/kparam_30/"
#     output:
#         outdir = "{indir}/clones/variants_init/knn/kparam_30/enriched_variants/donor{d}/params_results.tsv"
#     params:
#         donor = lambda wildcards: wildcards.d,
#         # Objective weights. order of the columns
#         weights = [1,0,0,1,-1, 1, 1], #weights=lambda wildcards: config["enriched_vars"]["params"]["objective_weights"][wildcards.objective_id]
#         objectives_l = ["variants_with_clone_norm_by_1_over_nclones_with_variant",
#                         "max_clone_ncells_over_nclones", "max_clone_ncells_over_ncells",
#                         "pct_thresh","other_pct_thresh", "n_vars", "obj_nclones_more_than_one_unique"],
#         ncpus=8, #config["ncpus"]
#         topn=16
#     notebook:
#         join(ROOT_DIR, "/workflow/notebooks/clone_vars/table_distinguishing_vars_multi_a_and_pctthresh.ipynb")
