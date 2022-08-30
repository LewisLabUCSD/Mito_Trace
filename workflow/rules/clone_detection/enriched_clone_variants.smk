from src.config import ROOT_DIR
from os.path import join, dirname
import numpy as np
import pandas as pd
import os
from icecream import ic
verbose = config.get("verbose", False)
if verbose:
    ic.enable()
else:
    ic.disable()



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
#ic('condition types for clonal shift', condition)

gff = config["genome_path"][config["genome"]]["gff"]
dendro_d = params_clones['params']['dendro_thresh']

objectives_l = params_clones["params"]["clone_barcodes_optimization"]["objectives"]
weights_cfg = params_clones["params"]["clone_barcodes_optimization"]["weights_cfg"]
constraint_ids = params_clones["params"]["clone_barcodes_optimization"]["constraints"]

#ic("weights_cfg", weights_cfg)
weights_ids = list(weights_cfg.keys())

#ic("weight ids", weights_ids)

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
        params_f = expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{{d}}/objs_{wt}/params.csv", wt = weights_ids)
    params:
        outdir = lambda wildcards, output: dirname(dirname(output.params_f[0])), #2 levels up
        indir = lambda wildcards, input: dirname(input.clone_cells),
        donor = lambda wildcards: wildcards.d,
        # Objective weights. order of the columns
        weights_cfg = weights_cfg,
        #weights = [1,1,1,1,-1, 1, 1], #weights=lambda wildcards: config["enriched_vars"]["params"]["objective_weights"][wildcards.objective_id]
        objectives_l = objectives_l, #["variants_with_clone_norm_by_1_over_nclones_with_variant", "max_clone_ncells_over_ncells","pct_thresh","other_pct_thresh", "n_vars", "obj_nclones_more_than_one_unique"],
        ncpus=12, #config["ncpus"]
        topn=8,
        to_test=False
    threads: 17
    priority: 8
    #script: join(ROOT_DIR, "workflow/notebooks/clone_vars/optimization_run.py.py")
    log:
        notebook = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/output.ipynb",
    notebook:
        join(ROOT_DIR, "workflow/notebooks/clone_vars/optimization_run_weights_wrapper.py.ipynb")


def get_constraints(wildcards):
    w = wildcards
    if w.cnstr == "None":
        return {"min_pct_thresh":0.2, "max_other_pct_thresh":0.8,
                "min_af_thresh": 0.005}
    elif w.cnstr == "min_pct_0.5":
        return {"min_pct_thresh":0.5, "max_other_pct_thresh":0.8,
                "min_af_thresh": 0.005}
    elif w.cnstr == "min_both_pct_0.5":
        return {"min_pct_thresh":0.5, "max_other_pct_thresh":0.5,
                "min_af_thresh": 0.005}
    else:
        raise ValueError("constraints property not set correctly")




# checkpoint run_variants_params_constraints:
#     """ Run the parameter optimization rule.
#
#         Get the best result and output the cells meta, af and dp.
#         For now also requires the labelled cells_meta file
#         for filtering later.
#     """
#     input:
#         clone_cells = "{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
#         anno_cells_meta_f = get_anno_cells #expand("{{outdir}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/gff_{gff}/annotation_clones/se_cells_meta_labels.tsv",
#                                    #gff=gff) # /cells_meta_labels.tsv"
#     output:
#         complete = expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{{d}}/objs_{wt}_constraints_{{cnstr}}/.complete", wt = weights_ids),
#         #params_f = expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{{d}}/objs_{wt}_constraints_{{cnstr}}/params.csv", wt = weights_ids),
#         #objs_norm_f = expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{{d}}/objs_{wt}_constraints_{{cnstr}}/objectives_norm.csv", wt = weights_ids)
#     params:
#         outdir = lambda wildcards, output: dirname(dirname(output.complete[0])), #2 levels up
#         indir = lambda wildcards, input: dirname(input.clone_cells),
#         donor = lambda wildcards: wildcards.d,
#         # Objective weights. order of the columns
#         weights_cfg = weights_cfg,
#         #weights = [1,1,1,1,-1, 1, 1], #weights=lambda wildcards: config["enriched_vars"]["params"]["objective_weights"][wildcards.objective_id]
#         objectives_l = objectives_l, #["variants_with_clone_norm_by_1_over_nclones_with_variant", "max_clone_ncells_over_ncells","pct_thresh","other_pct_thresh", "n_vars", "obj_nclones_more_than_one_unique"],
#         constraints = get_constraints,
#         constraint_name = lambda wildcards: f"constraints_{wildcards.cnstr}",
#         ncpus=12, #config["ncpus"]
#         topn=8,
#         to_test=False
#     threads: 17
#     priority: 8
#     #script: join(ROOT_DIR, "workflow/notebooks/clone_vars/optimization_run.py.py")
#     log:
#         notebook = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/{cnstr}_output.ipynb",
#     notebook:
#         join(ROOT_DIR, "workflow/notebooks/clone_vars/optimization_run_constraints_weights_wrapper.py.ipynb")
#
# rule optim_results:
#     input:
#         complete = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objs_{wt}_constraints_{cnstr}/.complete",
#         #params_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objs_{wt}_constraints_{cnstr}/params.csv",
#         #anno_cells_meta_f = get_anno_cells,
#         clone_cells = "{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
#     output:
#         #note = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/out.ipynb",
#         cells_meta_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objs_{wt}_constraints_{cnstr}/cells_meta.tsv",
#         best_params_f ="{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objs_{wt}_constraints_{cnstr}/best_params.csv",
#         best_params_clone_vars_f ="{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objs_{wt}_constraints_{cnstr}/best_params_filt_clone_vars.csv",
#         af_f ="{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objs_{wt}_constraints_{cnstr}/af.tsv",
#         dp_f ="{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objs_{wt}_constraints_{cnstr}/dp.tsv",
#         top_clustered = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objs_{wt}_constraints_{cnstr}/top_param_group_results.pdf",
#         top_table = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objs_{wt}_constraints_{cnstr}/top_param_groups_clone_vars.pdf",
#         top_clustered_table = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objs_{wt}_constraints_{cnstr}/top10perc_param_group_clusters.pdf",
#     priority: 10
#     params:
#         indir = lambda wildcards, input: dirname(input.complete),
#         outdir = lambda wildcards, output: dirname(output.top_clustered),
#         af_indirs = lambda wildcards, input: dirname(input.clone_cells),
#         #note = join(ROOT_DIR, "workflow/notebooks/clone_vars/optimization_results.ipynb"),
#         weights = lambda wildcards: weights_cfg[wildcards.wt], #[1,1,1,1,-1, 1, 1],
#         objectives_l = objectives_l, # ["variants_with_clone_norm_by_1_over_nclones_with_variant", "max_clone_ncells_over_ncells", "pct_thresh","other_pct_thresh", "n_vars", "obj_nclones_more_than_one_unique"],
#     log:
#         notebook = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objs_{wt}_constraints_{cnstr}/out.ipynb"
#     notebook:
#         join(ROOT_DIR, "workflow/notebooks/clone_vars/multi_cluster_optimization_results.py.ipynb")
#
# def get_multi_results(wildcards):
#     # note 1: ck_output is the same as OUTDIR, but this requires
#     # the checkpoint to complete before we can figure out what it is!
#
#     # note 2: checkpoints will have attributes for each of the checkpoint
#     # rules, accessible by name. Here we use make_some_files
#     ck_output = checkpoints.run_variants_params_constraints.get(**wildcards).output[0]
#     ic('ck_output', ck_output)
#     #SMP, = glob_wildcards(os.path.join(ck_output, "cells_meta.tsv"))
#     return glob_wildcards(os.path.join(dirname(ck_output), "cells_meta.tsv"))
#     #return expand(os.path.join(ck_output, "{SAMPLE}.txt"), SAMPLE=SMP)
#
#
# rule merge_donor_optim:
#     input:
#         cells_meta_f = get_multi_results,
#         #cells_meta_f = expand("{{outdir}}/enriched_barcodes/objs_{{wt}}_constraints_{{cnstr}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{d}/cells_meta.tsv", d=np.arange(config["N_DONORS"])),
#         # af_f =expand("{{outdir}}/enriched_barcodes/objs_{{wt}}_constraints_{{cnstr}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{d}/af.tsv", d=np.arange(config["N_DONORS"])),
#         # dp_f =expand("{{outdir}}/enriched_barcodes/objs_{{wt}}_constraints_{{cnstr}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{d}/dp.tsv", d=np.arange(config["N_DONORS"])),
#     output:
#         cells_meta_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
#         # af_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/af.tsv",
#         # dp_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/dp.tsv",
#     params:
#         indir = lambda wildcards, input: dirname(input.cells_meta_f),
#         outdir = lambda wildcards, output: dirname(output.cells_meta_f)
#     run:
#         pd.concat([pd.read_csv(join(x, "cells_meta.tsv"), sep="\t", index_col=0) for x in params["indir"]]).to_csv(output["cells_meta_f"], sep="\t")
#         pd.concat([pd.read_csv(join(x, "af.tsv"), sep="\t", index_col=0) for x in params["indir"]], axis=0).to_csv(join(params["outdir"], "af.tsv"), sep="\t")
#         pd.concat([pd.read_csv(join(x, "dp.tsv"), sep="\t", index_col=0) for x in params["indir"]], axis=0).to_csv(join(params["outdir"], "dp_f"), sep="\t")

rule run_variants_params_constraints:
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
        params_f = expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{{d}}/objs_{wt}_constraints_{{cnstr}}/params.csv", wt = weights_ids),
        objs_norm_f = expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{{d}}/objs_{wt}_constraints_{{cnstr}}/objectives_norm.csv", wt = weights_ids)
    params:
        outdir = lambda wildcards, output: dirname(dirname(output.params_f[0])), #2 levels up
        indir = lambda wildcards, input: dirname(input.clone_cells),
        donor = lambda wildcards: wildcards.d,
        # Objective weights. order of the columns
        weights_cfg = weights_cfg,
        #weights = [1,1,1,1,-1, 1, 1], #weights=lambda wildcards: config["enriched_vars"]["params"]["objective_weights"][wildcards.objective_id]
        objectives_l = objectives_l, #["variants_with_clone_norm_by_1_over_nclones_with_variant", "max_clone_ncells_over_ncells","pct_thresh","other_pct_thresh", "n_vars", "obj_nclones_more_than_one_unique"],
        constraints = get_constraints,
        constraint_name = lambda wildcards: f"constraints_{wildcards.cnstr}",
        ncpus=12, #config["ncpus"]
        topn=8,
        to_test=False
    threads: 17
    priority: 8
    #script: join(ROOT_DIR, "workflow/notebooks/clone_vars/optimization_run.py.py")
    log:
        notebook = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/{cnstr}_output.ipynb",
    notebook:
        join(ROOT_DIR, "workflow/notebooks/clone_vars/optimization_run_constraints_weights_wrapper.py.ipynb")


# checkpoint qc:
#     input:
#         objs_norm_f = expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{d}/objs_{{wt}}_constraints_{{cnstr}}/objectives_norm.csv",
#                              d=np.arange(config["N_DONORS"]))
#     output:
#         "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/qc_objs_{wt}_constraints_{cnstr}/qc.csv"
#     run:
#         cols = ["outdir", "variants", "kparam", "wt", "cnstr"]
#         vals = wildcards
#         ic('cols', cols)
#         ic('vals', vals)
#         obj_qc = pd.DataFrame(index = input["objs_norm_f"],
#                               columns=["objectives qc", "donor"] + cols)
#         for i, in_f in enumerate(list(input["objs_norm_f"])):
#             print('donor', i)
#             print('objs_norm head')
#             print(pd.read_csv(in_f).head())
#             if not (pd.read_csv(in_f).shape[0] < 1):
#                 obj_qc.loc[in_f, "objectives qc"] = True
#             else:
#                 obj_qc.loc[in_f, "objectives qc"] = False
#             obj_qc.loc[in_f, cols] = vals
#             obj_qc.loc[in_f, "donor"] = i
#         #ic(obj_qc.head())
#         obj_qc.index.name = "ID"
#         obj_qc.to_csv(output[0], index=True)
#         return


rule optim_results:
    input:
         #objs_norm_f = expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{{d}}/objs_{wt}_constraints_{{cnstr}}/objectives_norm.csv", wt = weights_ids)
        #params_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objs_{wt}_constraints_{cnstr}/params.csv",
        objs_norm_f = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/objs_{wt}_constraints_{cnstr}/objectives_norm.csv",
        anno_cells_meta_f = get_anno_cells,
        clone_cells = "{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
    output:
        #note = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/out.ipynb",
        cells_meta_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
        best_params_f ="{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/best_params.csv",
        best_params_clone_vars_f ="{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/best_params_filt_clone_vars.csv",
        af_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/af.tsv",
        dp_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/dp.tsv",
        top_clustered = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/top_param_group_results.pdf",
        top_table = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/top_param_groups_clone_vars.pdf",
        top_clustered_table = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/top10perc_param_group_clusters.pdf",
    priority: 10
    params:
        indir = lambda wildcards, input: dirname(input.objs_norm_f),
        outdir = lambda wildcards, output: dirname(output.top_clustered),
        af_indirs = lambda wildcards, input: dirname(input.clone_cells),
        #note = join(ROOT_DIR, "workflow/notebooks/clone_vars/optimization_results.ipynb"),
        weights = lambda wildcards: weights_cfg[wildcards.wt], #[1,1,1,1,-1, 1, 1],
        objectives_l = objectives_l, # ["variants_with_clone_norm_by_1_over_nclones_with_variant", "max_clone_ncells_over_ncells", "pct_thresh","other_pct_thresh", "n_vars", "obj_nclones_more_than_one_unique"],
    log:
        notebook = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/out.ipynb"
    notebook:
        join(ROOT_DIR, "workflow/notebooks/clone_vars/multi_cluster_optimization_results.py.ipynb")



rule tmp_optim:
    """ If no cells pass threshold for a donor, pass them all. Ideally
        should be removed from downstream tasksusing snakemake checkpoint,
        but some code issues arose. 
    
    """
    input:
        cells_meta_f = "{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
    output:
        cells_meta_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/temp/cells_meta.tsv",
        af_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/temp/af.tsv",
        dp_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/temp/dp.tsv",
    params:
        af_indir = lambda wildcards, input: join(dirname(input.cells_meta_f), "sc_af", f"donor{wildcards.d}"),
    run:
        print('wildcards tmp', wildcards["d"])
        df = pd.read_csv(input.cells_meta_f, sep="\t", index_col=0)
        df = df[df["donor"].astype(str) == str(wildcards["d"])]
        df["name"] = df.apply(lambda x: f'{x["donor"]}_{x["lineage"]}', axis=1)
        def inverse_string(x):
            spl = x.split("_")
            return f"{spl[1]}_{spl[0]}"
        df["ID"] = [inverse_string(x) for x in df.index]
        print('df shape')
        print(df.shape)
        df.to_csv(output.cells_meta_f, sep="\t")
        af_f = join(params.af_indir, "af.tsv")
        dp_f = join(params.af_indir, "dp.tsv")
        # transpose the af and dp DataFrames
        pd.read_csv(af_f, sep="\t", index_col=0).transpose().to_csv(output.af_f, sep="\t", index=True)
        pd.read_csv(dp_f, sep="\t", index_col=0).transpose().to_csv(output.dp_f, sep="\t", index=True)
        #outdir = dirname(output.cells_meta_f)
        #os.system(f"cp {af_f} {dp_f} {outdir}")



def check_donor(wildcards):
    w = wildcards
    print('w donor', w)
    in_f = f"{w.outdir}/enriched_barcodes/clones/variants_{w.variants}/knn/kparam_{w.kparam}/donor{w.d}/objs_{w.wt}_constraints_{w.cnstr}/objectives_norm.csv"
    print('in_f', in_f)
    print(pd.read_csv(in_f).shape[0] )
    if not (pd.read_csv(in_f).shape[0] < 1):
        in_d = f"{w.outdir}/enriched_barcodes/objs_{w.wt}_constraints_{w.cnstr}/clones/variants_{w.variants}/knn/kparam_{w.kparam}/donor{w.d}/"
    else:
        in_d = f"{w.outdir}/enriched_barcodes/objs_{w.wt}_constraints_{w.cnstr}/clones/variants_{w.variants}/knn/kparam_{w.kparam}/donor{w.d}/temp/"

    return [join(in_d,"cells_meta.tsv"), join(in_d,"af.tsv"), join(in_d,"dp.tsv")]
#        af_f ="{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/af.tsv",
        #dp_f ="{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/dp.tsv",



rule get_donor_cells:
    input:
        check_donor
    output:
        "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/single_donors/donor{d}/cells_meta.tsv",
        "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/single_donors/donor{d}/af.tsv",
        "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/single_donors/donor{d}/dp.tsv"
    params:
        outdir = lambda wildcards, output: dirname(output[0])
    shell:
        "cp {input} {params.outdir}/"


# def qc_filt(wildcards):
#     print('qc files', checkpoints.qc.get(**wildcards).output)
#     qc = pd.read_csv(checkpoints.qc.get(**wildcards).output[0], index_col=0)
#     qc = qc[qc["objectives qc"]]
#     return qc


# def get_results(wildcards):
#     w = wildcards
#     ic('w results')
#     ic(w)
#     outdir, wt, cnstr, variants, kparam = w
#     ic('outdir', outdir)
#     #qc = qc_filt(w)
#     #print('checkpoints no wiildcards', checkpoints.qc.get().output[0])
#     qc = pd.read_csv(checkpoints.qc.get(**wildcards).output[0], index_col=0)
#     print('qc get_results outdir', qc["outdir"].values)
#
#     qc = qc[qc["objectives qc"]]
#     #expand("{{outdir}}/enriched_barcodes/objs_{{wt}}_constraints_{{cnstr}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{d}/cells_meta.tsv", d=np.arange(config["N_DONORS"])),
#     files = []
#     print('qc', qc)
#     for ind, val in qc.iterrows():
#         curr_name = f"{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{val['donor']}/cells_meta.tsv"
#         files.append(curr_name)
#     return expand("{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{donor}/cells_meta.tsv",
#                   zip, outdir=qc["outdir"].values, wt=qc["wt"].values, cnstr=qc["cnstr"].values, variants=qc["variants"].values,
#                   kparam=qc["kparam"].values, donor=qc["donor"].values)



rule merge_donor_optim:
    input:
        #cells_meta_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/single_donors/donor{d}/cells_meta.tsv"
        #cells_meta_f = get_results
        cells_meta_f = expand("{{outdir}}/enriched_barcodes/objs_{{wt}}_constraints_{{cnstr}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/single_donors/donor{d}/cells_meta.tsv", d=np.arange(config["N_DONORS"])),
        af_f = expand("{{outdir}}/enriched_barcodes/objs_{{wt}}_constraints_{{cnstr}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/single_donors/donor{d}/af.tsv", d=np.arange(config["N_DONORS"])),
        dp_f = expand("{{outdir}}/enriched_barcodes/objs_{{wt}}_constraints_{{cnstr}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/single_donors/donor{d}/dp.tsv", d=np.arange(config["N_DONORS"])),
        # af_f =expand("{{outdir}}/enriched_barcodes/objs_{{wt}}_constraints_{{cnstr}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{d}/af.tsv", d=np.arange(config["N_DONORS"])),
        # dp_f =expand("{{outdir}}/enriched_barcodes/objs_{{wt}}_constraints_{{cnstr}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/donor{d}/dp.tsv", d=np.arange(config["N_DONORS"])),
    output:
        cells_meta_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
        af_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/af.tsv",
        # dp_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/dp.tsv",
    params:
        #indir = lambda wildcards, input: dirname(input.cells_meta_f),
        outdir = lambda wildcards, output: dirname(output.cells_meta_f)
    #priority: 10
    run:
        print("input merge")
        print(input.cells_meta_f)
        print(type(input["cells_meta_f"]))
        print('cells meta', input["cells_meta_f"])
        pd.concat([pd.read_csv(join(dirname(x), "cells_meta.tsv"), sep="\t", index_col=0) for x in input["cells_meta_f"]], axis=0).to_csv(output["cells_meta_f"], sep="\t")
        pd.concat([pd.read_csv(join(dirname(x), "af.tsv"), sep="\t", index_col=0) for x in input["cells_meta_f"]], axis=0).fillna(0).to_csv(join(params["outdir"], "af.tsv"), sep="\t")
        pd.concat([pd.read_csv(join(dirname(x), "dp.tsv"), sep="\t", index_col=0) for x in input["cells_meta_f"]], axis=0).fillna(0).to_csv(join(params["outdir"], "dp.tsv"), sep="\t")

        # pd.concat([pd.read_csv(x, sep="\t", index_col=0) for x in input.cells_meta_f]).to_csv(output.cells_meta_f, sep="\t")
        # pd.concat([pd.read_csv(x, sep="\t", index_col=0) for x in input.af_f], axis=0).to_csv(output.af_f, sep="\t")
        # pd.concat([pd.read_csv(x, sep="\t", index_col=0) for x in input.dp_f], axis=0).to_csv(output.dp_f, sep="\t")



rule barcodes_btwnClones_dendro:
    input:
        cells_meta = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
        #af_note = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/af.tsv",
    output:
        note = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/btwnClones_dendro/donor{d}.ipynb",
        mean = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/btwnClones_dendro/donor{d}.mean.csv",
        res = report(multiext("{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/btwnClones_dendro/donor{d}.",
                                  "clones_dendro.csv", "dendrogram_pvals.txt",
                                  "dendro.NoCondition.max2.AF.png"),
                         category="lineage")
    params:
        note = join(ROOT_DIR, "workflow", "notebooks", "clone_af_dendrograms", "MT_btwnClones_enrichedBarcodes_dendro_dynamic.ipynb"),
        indir = lambda wildcards, input: dirname(input.cells_meta),
        outdir = lambda wildcards, output: dirname(output.note),
        #dendro_thresh = 0.6
    shell: "papermill -p INDIR {params.indir} -p OUTDIR {params.outdir} -p DONOR {wildcards.d} {params.note} {output.note}"

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


name_map = {"clones":"name", "dendro_bc": "den_clust", "mt_as_clones": "Variants", "mt_as_clones_dendro": "den_clust"}

module cloneShiftMod:
    snakefile: "../clonal_shift.smk"
    config: config

#use rule * from cloneShiftMod as cloneshifttwo_*

module reprClonesMod:
    snakefile: "../repr_clones.smk"
    config: config

#use rule * from reprClonesMod as reprcl_*


rule counts_clones:
    input:
        #se_meta = "{outdir}/annotation_clones/se_cells_meta.tsv",
        se_meta = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/single_donors/donor{d}/cells_meta.tsv" #"{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
        #se_meta = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
    output:
        out = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/annotation_clones/clone_counts/donor{d}/clone_counts.barplot_conditions.png",
        note = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/annotation_clones/clone_counts/donor{d}/counts_clones.ipynb"
    params:
        outdir = lambda wildcards, output: dirname(output.note),
        note = join(ROOT_DIR, "workflow/notebooks/lineage_clones/python_clone_cell_counts.ipynb"),
        sample_names = ",".join(config["samples"].index),
        min_cell = config["annotation_clones"]["params"]["clone_sizes_min_cell"]
    shell:
        "papermill -p se_cells_meta_f {input.se_meta} -p outdir {params.outdir} -p sample_names {params.sample_names} -p min_cell {params.min_cell} {params.note} {output.note}"

#cells_meta_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",

## Clonal shift rules:
## a) rule clonalshift_clones  and clonalshift_dendro_bc
## b) rule finalize
#######################################################################
use rule clonalshift_clones from cloneShiftMod as cloneshifttwo_clonalshift_clones with:
    input:
        se_meta = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv", #"{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
        # se_meta = expand("{{outdir}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/annotation_clones/se_cells_meta_labels.tsv",
        #              gff=config["gff"])
    output:
        note = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clonal_shifts/variants_{variants}/clones/results/{condition}/knn/kparam_{kparam}/output.ipynb"

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
        noInput = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clonal_shifts/variants_{variants}/clones/results/noInput/knn/kparam_{kparam}/output.ipynb",
        input = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clonal_shifts/variants_{variants}/clones/results/inputOnly/knn/kparam_{kparam}/output.ipynb"
    output:
        note = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clonal_shifts/variants_{variants}/donors/donor{d}/clones/knn_kparam_{kparam}/output_compare.ipynb"
    params:
        clone_col = name_map["clones"], #lambda wildcards: name_map[wildcards.cloneShift_method]
        noInput_indir= lambda wildcards, input: dirname(input.noInput),
        input_indir = lambda wildcards, input: dirname(input.input),
        outdir = lambda wildcards, output: dirname(output.note),
        donor = lambda wildcards: wildcards.d,
        script = join(ROOT_DIR, "workflow/notebooks/clonal_shifts/combine_conditions_hypergeometric.ipynb"), #get_script,



use rule single_input_and_culture_cl from cloneShiftMod as cloneshifttwo_single_input_and_culture_cl with:
    input:
        cond = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clonal_shifts/variants_{variants}/clones/results/{condition}/knn/kparam_{kparam}/output.ipynb",
    output:
        note = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clonal_shifts/variants_{variants}/donors/donor{d}/clones/knn_kparam_{kparam}/output_single_{condition}_compare.ipynb"
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
        return f"{w.outdir}/enriched_barcodes/objs_{w.wt}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/clones/knn_kparam_{w.kparam}/output_compare.ipynb"
    else:
        if "inputOnly" in condition:
            return f"{w.outdir}/enriched_barcodes/objs_{w.wt}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/clones/knn_kparam_{w.kparam}/output_single_inputOnly_compare.ipynb"
        elif "noInput" in condition:
            return f"{w.outdir}/enriched_barcodes/objs_{w.wt}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/clones/knn_kparam_{w.kparam}/output_single_noInput_compare.ipynb"
    return


def get_compare_output_constraints(wildcards):
    w = wildcards
    if len(condition) > 1:
        return f"{w.outdir}/enriched_barcodes/objs_{w.wt}_constraints_{w.cnstr}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/clones/knn_kparam_{w.kparam}/output_compare.ipynb"
    else:
        if "inputOnly" in condition:
            return f"{w.outdir}/enriched_barcodes/objs_{w.wt}_constraints_{w.cnstr}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/clones/knn_kparam_{w.kparam}/output_single_inputOnly_compare.ipynb"
        elif "noInput" in condition:
            return f"{w.outdir}/enriched_barcodes/objs_{w.wt}_constraints_{w.cnstr}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/clones/knn_kparam_{w.kparam}/output_single_noInput_compare.ipynb"
    return


use rule mv_tmp_compare from cloneShiftMod as cloneshifttwo_mv_tmp_compare with:
    input:
        get_compare_output_constraints
    output:
        note = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clonal_shifts/variants_{variants}/donors/donor{d}/clones/knn_kparam_{kparam}/output.ipynb"

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
         expand("{{outdir}}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clonal_shifts/enriched_barcodes/variants_{{variants}}/dendro_bc/results/{condition}/knn/kparam_{kparam}/output.ipynb",
                condition=condition, kparam=params_clones["knn"]["params"]["resolution"], wt = weights_ids, cnstr=constraint_ids),
         expand("{{outdir}}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clonal_shifts/variants_{{variants}}/clones/results/{condition}/knn/kparam_{kparam}/output.ipynb",
                condition=condition, kparam=params_clones["knn"]["params"]["resolution"], wt = weights_ids, cnstr=constraint_ids),
    output:
          "{outdir}/enriched_barcodes/clonal_shifts/variants_{variants}/.finalize_complete.txt"
    #shell: "touch {output}"



## Indclones
## a) rule get_clone_cells and get_clone_dendro_bc
## b) top_clone_complete
#######################################################################
use rule get_clone_cells from reprClonesMod as reprcl_get_clone_cells with:
    input:
        #se_meta = "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/se_cells_meta_labels.tsv"
        se_cells_meta_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv", #"{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
        # se_meta = expand("{{outdir}}/enriched_barcodes/clones/variants_{{variants}}/knn/kparam_{{kparam}}/annotation_clones/cells_meta.tsv",
        #              gff=config["gff"])
    output:
        cells_meta = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv"



use rule get_repr_clones_cl from reprClonesMod as reprcl_get_repr_clones_cl with:
    input:
        hyper = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clonal_shifts/variants_{variants}/donors/donor{d}/clones/knn_kparam_{kparam}/output.ipynb", #get_hypergeo,
        cells="{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv"
    output:
        note = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/clones_ranked/output_rank_ncells.ipynb",
        clone_order = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/clones_ranked/representative_cloneID.txt"



use rule top_clone_mt_variants_cl from reprClonesMod as reprcl_top_clone_mt_variants_cl with:
    input:
        #se =  expand("{{outdir}}/anno_multiplex/gff_{gff}/SE.rds", gff=config["gff"]),
        cells = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv",
        clone_order_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/clones_ranked/representative_cloneID.txt",
        mt_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/af.tsv",
    output:
        note = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/{clone_order}_clone_mt_variants.ipynb",


# use rule top_clone_mt_variants_mt from reprClonesMod as reprcl_clone_mt_variants_mt with:
#     input:
#         se =  expand("{{outdir}}/anno_multiplex/gff_{gff}/SE.rds", gff=config["gff"]),
#         cells = "{outdir}/enriched_barcodes/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv",
#         mt_f =  "{outdir}/enriched_barcodes/clones/variants_{variants}/knn/kparam_{kparam}/af.tsv",
#     output:
#         note = "{outdir}/enriched_barcodes/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_clones/top/{clone_order}_clone_mt_variants.ipynb",
#

use rule top_clone_umap_overlay from reprClonesMod as reprcl_top_clone_umap_overlay with:
    input:
        se =  expand("{{outdir}}/anno_multiplex/gff_{gff}/SE.rds", gff=config["gff"]),
        cells = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv",
        clone_order_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/clones_ranked/representative_cloneID.txt",
    output:
        note = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/top_umap_overlay.ipynb"
    params:
        script = join(ROOT_DIR, "workflow/notebooks/individual_clones/top_individual_clones_umap_overlap.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),


def get_top_clone_hypergeoSig_script(wildcards):
    return join(ROOT_DIR, "workflow/notebooks/individual_clones/oneCondition_top_individual_clone_lineage_Sig_hypergeo.ipynb")

def get_param_top_clone_hypergeoSig_condition():
    if len(condition) == 1:
        return condition[0]
    else:
        return "None"


use rule top_clone_hypergeoSig_cl from reprClonesMod as reprcl_top_clone_hypergeoSig_cl with:
    input:
        indir = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clonal_shifts/variants_{variants}/donors/donor{d}/clones/knn_kparam_{kparam}/output.ipynb", #get_hypergeo,
        cloneIDs_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/clones_ranked/representative_cloneID.txt",
        #cloneID = "{outdir}/enriched_barcodes/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/cloneIDs/{cloneID}.txt"
    output:
        note = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/top/top_hypergeo_sig.ipynb",


use rule top_clone_lineage_count from reprClonesMod as reprcl_top_clone_lineage_count with:
    input:
        cells = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv",
        clone_order_f = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/clones_ranked/representative_cloneID.txt",
    output:
        note = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/top_clone_lineage_count.ipynb",


use rule top_merge from reprClonesMod as reprcl_top_merge with:
    input:
        ins = multiext("{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/","top_umap_overlay.ipynb", "top_hypergeo_sig.ipynb", "top_clone_lineage_count.ipynb", "top_clone_mt_variants.ipynb"),
    output:
        out_f = multiext("{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/","clone_shift_combine.svg"),
        note = "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/merge_panels.ipynb"




# def get_repr_results(wildcards, col_names=("outdir", "wt", "cnstr", "variants", "kparam", "donor")):
#     w = wildcards
#     f=["top_umap_overlay.ipynb", "top_hypergeo_sig.ipynb", "top_clone_lineage_count.ipynb", "all_clone_mt_variants.ipynb", "clone_shift_combine.svg"],
#     clone_subset=[]
#     ic('w repr results')
#     ic(w)
#     outdir, wt, cnstr, clone_subset = w
#
#     ic('outdir', outdir)
#     qc = pd.read_csv(checkpoints.qc.get(**wildcards).output[0], index_col=0)
#     qc = qc[qc["objectives qc"]]
#     ic(qc.head())
#     #qc = qc_filt(w)
#     files = []
#     for ind, val in qc.iterrows():
#         outdir, wt, cnstr, variants, kparam, d = val[col_names].values
#         curr_name = f"{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/variants_{val['variants']}/knn/kparam_{val['kparam']}/donor{val['donor']}/cells_meta.tsv",
#         files.append(curr_name)
#     ic(files)
#     return expand("{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/variants_{variants}/knn/kparam_{kparam}/donor{donor}/cells_meta.tsv",
#                   zip, outdir=qc["outdir"].values, wt=qc["wt"].values, cnstr=qc["cnstr"].values, variants=qc["variants"].values,
#                   kparam=qc["kparam"].values, donor=qc["donor"].values)
#
#     #return files

# def get_repr_clone_out(wildcards):
#     w = wildcards
#     all_files = []
#     for variants in params_clones["variants"]:
#         for kparam in params_clones["knn"]["params"]["resolution"]:
#             for d in np.arange(config["N_DONORS"]):
#                 for f in ["top_umap_overlay.ipynb", "top_hypergeo_sig.ipynb", "top_clone_lineage_count.ipynb", "all_clone_mt_variants.ipynb", "clone_shift_combine.svg"]:
#                     in_f = ""
#                     all_files.append(f"{w.outdir}/enriched_barcodes/objs_{w.wt}_constraints_{w.cnstr}/{w.clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/{f}")
#     return


use rule top_clone_complete from reprClonesMod as reprcl_top_clone_complete with:
    input:
        expand("{{outdir}}/enriched_barcodes/objs_{{wt}}_constraints_{{cnstr}}/{{clone_subset}}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/top/{f}",
              f=["top_umap_overlay.ipynb", "top_hypergeo_sig.ipynb", "top_clone_lineage_count.ipynb", "all_clone_mt_variants.ipynb", "clone_shift_combine.svg"],
              d=np.arange(config["N_DONORS"]), variants=[x for x in params_clones["variants"] if x != "simple"],
              kparam=params_clones["knn"]["params"]["resolution"],)
        #get_repr_results
    output:
        "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/.top_clones.txt"



def get_obj_names():
    names = ""
    for wt in weights_ids:
        for cnstr in constraint_ids:
            if names != "":
                names = names + "..."
            names = names + f"objs_{wt}_constraints_{cnstr}"
    return names

def get_opt_top(wildcards):
    w = wildcards
    return f"single_clones/donor{w.d}/cloneMethod_variants_{w.variants}_knn_resolution_{w.kparam}/clonalShift_method_clones/top/",


rule compare_clone_var_tables:
    input:
        cells_meta_f = expand("{{outdir}}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/cells_meta.tsv",
                              wt=weights_ids, cnstr=constraint_ids),
    output:
        out_f = "{outdir}/clones/comparisons/variants_{variants}/knn/kparam_{kparam}/donor{d}/clones_dendro.svg",
        note = "{outdir}/clones/comparisons/variants_{variants}/knn/kparam_{kparam}/donor{d}/compare_clone_tables.ipynb",
        #variants_dirs = ["variants_init", "variants_simpleUnion", "variants_preFilterImpute"]
    params:
        donor = lambda wildcards: wildcards.d,
        indir = lambda wildcards: wildcards.outdir,
        outdir = lambda wildcards, output: dirname(output.out_f),
        prm_variants = lambda wildcards: f"variants_{wildcards.variants}",
        prm_clone_k = lambda wildcards: f"kparam_{wildcards.kparam}",
        prms_enriched_barcode = get_obj_names(),
        prm_opt_top = get_opt_top,
        gff_id = f"gff_{gff}",
        note = join(ROOT_DIR, "workflow/notebooks/aggregate_experiments/concat_single_donor_figures.py.ipynb")
    shell: "papermill -p indir {params.indir} -p outdir {params.outdir} -p donor {params.donor} -p gff_id {params.gff_id} -p prm_variants {params.prm_variants} -p prm_clone_k {params.prm_clone_k} -p prms_enriched_barcode {params.prms_enriched_barcode} -p  prm_opt_top {params.prm_opt_top} -p gff_id {params.gff_id} {params.note} {output.note}"

# rule finalize:
#     input:
#         expand("{{outdir}}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/btwnClones_dendro/donor{d}.ipynb",
#                 variants=[x for x in params_clones["variants"] if x != "simple"],
#                 kparam=params_clones["knn"]["params"]["resolution"],
#                 wt = weights_ids, d=np.arange(config["N_DONORS"])),
#         expand("{{outdir}}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/.top_clones.txt", wt = weights_ids, clone_subset=["repr_clones", "single_clones"]),
#         expand("{{outdir}}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/annotation_clones/clone_counts/donor{d}/counts_clones.ipynb",
#                variants=[x for x in params_clones["variants"] if x != "simple"],
#                kparam=params_clones["knn"]["params"]["resolution"],
#                wt = weights_ids, d=np.arange(config["N_DONORS"]))
#     output: final = "{outdir}/enriched_barcodes/.finalize"
#     shell: "touch {output}"


############################################
## Compare methods:
############################################
rule compare_clonevar_methods_get_params:
    input:
        expand("{{outdir}}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
               variants=[x for x in params_clones["variants"] if x != "simple"], kparam=params_clones["knn"]["params"]["resolution"],
               wt=weights_ids, cnstr=constraint_ids),
    output:
        meta_params = "{outdir}/enriched_barcodes/comparisons/params_meta.csv",
    # params:
    #     pre_outdir = lambda wildcards: wildcards.outdir
    run:
        import itertools
        df = pd.DataFrame(index=input, columns=["variants", "kparam", "objective", "constraints"])
        params_iter = list(itertools.product(params_clones["variants"],params_clones["knn"]["params"]["resolution"],
                                              weights_ids, constraint_ids))
        ic('params_iter')
        ic(params_iter[:2])
        for x in params_iter:
            variants, kparam, wt, cnstr = x
            df.loc[f"{wildcards['outdir']}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv"] = [variants, kparam, wt, cnstr]
        ic(df.head())
        df.to_csv(output["meta_params"])


rule compare_clonevar_methods:
    input:
        meta_params = "{outdir}/enriched_barcodes/comparisons/params_meta.csv",
    output:
        note = "{outdir}/enriched_barcodes/comparisons/comparisons.ipynb",
        #fig = report(multiext("{outdir}/enriched_barcodes/comparisons/", "methods_nAgreeNorm.png", "methods_nAgreeNorm_agg.png", "methods_n_T_T_agg.png", "methods_nTTNorm_agg.png"), category="Methods"),
    params:
        note = join("src", "clones_compare", "distance_matrix_enriched.ipynb"),
        outdir = lambda wildcards, output: dirname(output[0]),
        # prefix = config["prefix"],
        # cfg_outdir = config["outdir"],
        # params_f = config["config"],
        # all = lambda wildcards, input: ",".join(x for x in input),
    shell:
        "papermill -p meta_params_f {input.meta_params} -p outdir {params.outdir} {params.note} {output.note}"
        #"papermill -p all_files {params.all} -p outdir {params.outdir} -p cfg_outdir {params.cfg_outdir} -p prefix {params.prefix} -p params_f {params.params_f} {params.note} {output.note}"





# def get_col_values(w, col_names):
#     col_values = []
#     for ind, val in enumerate(w):
#         col_values.append(col_names[ind])
#     return col_values

#
# def get_dendro_filt(wildcards, col_names):
#     w = wildcards
#     qc = qc_filt(w)
#     col_values = get_col_values(w, col_names)
#     outdir, wt, cnstr, variants, kparam = w
#     files = []
#     for ind, val in qc.iterrows():
#         curr_name = f"{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{val['donor']}/cells_meta.tsv",
#         files.append(curr_name)
#     return files
#
# def get_counts_clones(wildcards):
#     w = wildcards
#     qc = qc_filt(w)
#     return
#
#
# def get_filt_compare_clones(wildcards):
#     w = wildcards
#     qc = qc_filt(w)
#     #outdir, wt, cnstr, variants, kparam = w
#     files = []
#     for ind, val in qc.iterrows():
#         curr_name = f"{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/donor{val['donor']}/cells_meta.tsv",
#         files.append(curr_name)
#     return files




def param_names(qc, col_names=("outdir", "wt", "cnstr", "variants", "kparam", "donor")):
    ic('param_names')
    files = []
    for ind, val in qc.iterrows():
        outdir, wt, cnstr, variants, kparam, d = val[list(col_names)].values
        ic('outdir', outdir)
        curr_objs_dir = f"{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}"
        curr_clones_dir = join(curr_objs_dir, f"clones/variants_{variants}/knn/kparam_{kparam}")
        #clone_sub_dir = join(curr_objs_dir, clone_subset, ".top_clones.txt")
        clcounts_f = join(curr_clones_dir, f"annotation_clones/clone_counts/donor{d}/counts_clones.ipynb")
        dendro_f = join(curr_clones_dir, "barcodes", f"btwnClones_dendro/donor{d}.ipynb")
        clone_tables_f = join(outdir, f"clones/comparisons/variants_{variants}/knn/kparam_{kparam}/donor{d}/compare_clone_tables.ipynb")
        files.append(dendro_f)
        files.append(clcounts_f)
        files.append(clone_tables_f)
        #files.extend([dendro_f, clcounts_f, clone_tables_f])
    print("files summary", files)
    return files

#
# def summarize_donors(wildcards):
#     #qc = qc_filt(wildcards)
#     ic('wildcards summarize')
#     ic(wildcards)
#     ic('qc')
#     qc = pd.read_csv(checkpoints.qc.get(**wildcards).output[0], index_col=0)
#     qc = qc[qc["objectives qc"]]
#     ic('qc', qc.head())
#     return param_names(qc)


# rule finalize_donors:
#     input:
#         summarize_donors
#     output:
#         "{outdir}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/.summarize_donors",
#                # variants=[x for x in params_clones["variants"] if x != "simple"],
#                # kparam=params_clones["knn"]["params"]["resolution"],
#                # wt = weights_ids, cnstr = constraint_ids,

rule finalize_constraints:
    input:
        #summarize_donors,
        # expand("{{outdir}}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/.summarize_donors",
        #        variants=[x for x in params_clones["variants"] if x != "simple"],
        #        kparam=params_clones["knn"]["params"]["resolution"],
        #        wt = weights_ids, cnstr = constraint_ids),
        expand("{{outdir}}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/{clone_subset}/.top_clones.txt",
               wt = weights_ids, cnstr = constraint_ids, clone_subset=["repr_clones", "single_clones"]),

        expand("{{outdir}}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/btwnClones_dendro/donor{d}.ipynb",
               variants=[x for x in params_clones["variants"] if x != "simple"],
               kparam=params_clones["knn"]["params"]["resolution"],
               wt = weights_ids, cnstr = constraint_ids,
               d=np.arange(config["N_DONORS"])),
        expand("{{outdir}}/enriched_barcodes/objs_{wt}_constraints_{cnstr}/clones/variants_{variants}/knn/kparam_{kparam}/annotation_clones/clone_counts/donor{d}/counts_clones.ipynb",
               variants=[x for x in params_clones["variants"] if x != "simple"],
               kparam=params_clones["knn"]["params"]["resolution"],
               wt = weights_ids, cnstr = constraint_ids,
               d=np.arange(config["N_DONORS"])),
        expand("{{outdir}}/clones/comparisons/variants_{variants}/knn/kparam_{kparam}/donor{d}/compare_clone_tables.ipynb",
               variants=[x for x in params_clones["variants"] if x != "simple"],
               kparam=params_clones["knn"]["params"]["resolution"],
               d=np.arange(config["N_DONORS"])),

        "{outdir}/enriched_barcodes/comparisons/comparisons.ipynb",

    output: final = "{outdir}/enriched_barcodes/.finalize"
    shell: "touch {output}"



