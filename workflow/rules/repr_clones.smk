import pandas as pd
from src.config import ROOT_DIR
from src.utils.parse_config import read_config_file
import os
import numpy as np
from os.path import join, dirname
from snakemake.utils import min_version
min_version("6.0")

from icecream import ic
verbose = config.get("verbose", False)
if verbose:
    ic.enable()
else:
    ic.disable()

knn_clone_shift= ["clones", "dendro_bc"]
#mt_method = "mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}"
mt_clone_shift = ["mt_as_clones_dendro"]#, "mt_as_clones"]
best_p = config["mt_as_clones"]["best_params"]
params_clones = config["clones"]


dendro_d = config['clones']['params']['dendro_thresh']


sample_ids = config["samples"].index
condition = []
if "Input" in sample_ids:
    condition.append("inputOnly")
    if len(sample_ids) > 1:
        condition.append("noInput")
elif len(sample_ids) > 0:
    condition.append("noInput")



def load_cells_meta(cells_meta_f):
    ic(cells_meta_f)
    ic('here')
    try:
        cells_meta = pd.read_csv(cells_meta_f, sep="\t", index_col=0)
    except:
        cells_meta = pd.read_csv(cells_meta_f[0], sep="\t", index_col=0)
    #cells_meta = pd.read_csv(cells_meta_f, sep="\t", index_col=0)
    cells_meta = cells_meta.loc[~(cells_meta["name"]=="None")]
    if not "cluster_labels" in cells_meta.columns.values:
        cells_meta["clusterID"] = cells_meta["seurat_clusters"]
    else:
        cells_meta["clusterID"] = cells_meta["cluster_labels"]
    ic('cells_meta')
    ic(cells_meta.head())
    return cells_meta




""" Output line-separated cell IDs
Cell ID, Clone, Lineage, Condition
"""
rule get_clone_cells:
    input:
        se_cells_meta_f = expand("{{outdir}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/gff_{gff}/annotation_clones/se_cells_meta_labels.tsv", gff=config["gff"])
    output:
        cells_meta = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv"
    run:
        cells_meta = load_cells_meta(input.se_cells_meta_f)
        ic('cells_meta', cells_meta.head())
        cells_meta = cells_meta.rename({"name":"cloneID"}, axis=1)
        cells_meta = cells_meta.loc[cells_meta["donor"].astype(int) == int(wildcards.d)]
        ic('cells_meta', cells_meta.head())
        cells_meta[["cloneID", "clusterID","condition", "donor"]].to_csv(output.cells_meta, sep="\t")


rule get_clone_dendro_cells:
    input:
        barcodes_dir = expand("{{outdir}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/barcodes/btwnClones_dendro_dt_{dt}/donor{{d}}.clones_dendro.csv", dt=dendro_d),
        se_cells_meta_f = expand("{{outdir}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/gff_{gff}/annotation_clones/se_cells_meta_labels.tsv", gff=config["gff"])#"{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/se_cells_meta_labels.tsv"
    output:
        cells = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_dendro_bc/cells_meta.tsv"
    #output: "{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/{clone}/cells.txt"
    run:
        clone_col = "den_clust"
        cells_meta = load_cells_meta(input.se_cells_meta_f)
        #barcodes_in = pd.read_csv(join(input.barcodes_dir[0],f"donor{wildcards.d}.clones_dendro.csv"), index_col=0)
        barcodes_in = pd.read_csv(list(input.barcodes_dir)[0], index_col=0)
        ic('d', wildcards.d)
        #ic('barcodes_in')
        ic('a barcodes',barcodes_in.head())
        barcodes_in[clone_col] = str(wildcards.d) + "_" + barcodes_in[clone_col]
        ic('b',barcodes_in.head())
        cells_meta = cells_meta[cells_meta["name"].isin(barcodes_in.index)]
        cells_meta[clone_col] = cells_meta.apply(lambda x: barcodes_in.loc[x["name"], clone_col] , axis=1)
        ic('c', cells_meta.head())
        cells_meta = cells_meta.loc[cells_meta["donor"].astype(str) == str(wildcards.d)]

        cells_meta = cells_meta.rename({clone_col: "cloneID"}, axis=1)
        cells_meta[["cloneID", "clusterID","condition", "donor"]].to_csv(output.cells, sep="\t")


rule get_cells_mt_as_clones:
    input:
        clones = "{outdir}/mt_as_clones/variants_{variants}/bestparams/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/best_params_save.ipynb",
        se_cells_meta_f = expand("{{outdir}}/anno_multiplex/gff_{gff}/se_cells_meta_labels.tsv", gff=config["gff"])
    output:
        cells="{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_mt_as_clones/cells_meta.tsv"
        #"{outdir}/{clone_subset}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/{clone}/cells.txt"
    params:
        indir = lambda wildcards, input: dirname(input.clones),
        pct_thresh=0.8
    run:
        indir = params.indir
        best_params = pd.read_csv(join(indir,  "best_params.csv"))
        bin_d = pd.read_csv(join(indir,  f"donor_{wildcards.d}_binary_na.csv"), index_col=0)
        bin_d.columns = [f"cloneID_{x}" for x in bin_d.columns]
        labels_df = pd.read_csv(list(input.se_cells_meta_f)[0],sep="\t").set_index("ID")
        ic(labels_df.shape)
        labels_df = labels_df[~(labels_df["donor"]=='None')]
        ic(labels_df.shape)
        cells_meta = pd.merge(labels_df[["donor", "cluster_labels", "condition"]], bin_d, left_index=True, right_index=True, how="inner")
        cells_meta.index = [f'{x.split("_")[1]}_{x.split("_")[0]}'  for x in cells_meta.index]
        cells_meta = cells_meta.rename({"cluster_labels":"clusterID"}, axis=1)

        # Drop donor-specific variants (i.e. pct_thresh of cells have the clone)
        n_pct = params.pct_thresh*cells_meta.shape[0]
        cell_clones = cells_meta.drop(["donor", "clusterID", "condition"], axis=1).copy()
        donor_vars = cell_clones.loc[:,((cell_clones>0).sum()>n_pct)].columns.values
        #donor_vars = get_high_variants(cells_meta.drop(["donor", "clusterID", "condition"], axis=1), thresh=0, pct_thresh=0.8)
        cells_meta = cells_meta.loc[:, ~(cells_meta.columns.isin(donor_vars))]
        #donor_vars = get_high_variants(cells_meta, thresh=0.8, pct_thresh=0.9)
        ic(f"number of donor vars: {len(donor_vars)}")
        cells_meta.to_csv(output.cells, sep="\t")

# convert back to raw cell count from log2
#lin_mt_ncells["count"] = (np.ceil((2**lin_mt_ncells["count"])-1)).astype(int)

rule get_cells_mt_as_clones_dendro:
    input:
        clones = "{outdir}/mt_as_clones/variants_{variants}/dendro/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/dendro_mt_clones.ipynb",
        se_cells_meta_f = expand("{{outdir}}/anno_multiplex/gff_{gff}/se_cells_meta_labels.tsv", gff=config["gff"])
    output:
        cells="{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_mt_as_clones_dendro/cells_meta.tsv"
    params:
        indir = lambda wildcards, input: dirname(input.clones)
    run:
        indir = params.indir
        clone_col="den_clust"
        atac_col = "cluster_labels"
        labels_df = pd.read_csv(list(input.se_cells_meta_f)[0],sep="\t").reset_index().set_index("ID")
        ic(labels_df.shape)
        labels_df = labels_df[~(labels_df["donor"]=='None')]
        ic(labels_df.shape)
        labels_df.head()
        den_d = pd.read_csv(join(indir, f"don_{wildcards.d}_mt_dendro_clust.csv"), index_col=0)
        labels = labels_df[labels_df["donor"] == wildcards.d]
        cells_meta = pd.merge(left=den_d, right=labels[[atac_col]], left_index=True, right_index=True, how="inner")
        cells_meta.index = [f'{x.split("_")[1]}_{x.split("_")[0]}'  for x in cells_meta.index]
        cells_meta = cells_meta.rename({"Donor":"donor", "cluster_labels":"clusterID", clone_col:"cloneID"},axis=1)
        cells_meta.to_csv(output.cells, sep="\t")


rule get_cloneIDs:
    input:
        cells = "{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/cells_meta.tsv"
    output:
        clone_f = "{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/cloneIDs.txt",
    run:
        df = pd.read_csv(input.cells, sep="\t")
        if wildcards.clone_shift_method == "mt_as_clones":
            cloneIDs = [x for x in df.columns if "cloneID" in x ]
        else:
            cloneIDs = list(set(df["cloneID"]))
        ic('cloneIDs', cloneIDs)

        cloneIDs_str = "\n".join(cloneIDs)
        with open(output.clone_f, "w") as f:
            f.write(cloneIDs_str)


def get_hypergeo(wildcards):
    w = wildcards
    clshift = w.clone_shift_method
    indir = f"{w.outdir}/clonal_shifts/variants_{w.variants}/donors/donor{w.d}/{w.clone_shift_method}"
    #output.ipynb
    if  clshift == "dendro_bc" or clshift == "clones":
        return f"{indir}/knn_kparam_{w.kparam}/output.ipynb"
    elif clshift == "mt_as_clones_dendro" or clshift == "mt_as_clones":
        return f"{indir}/af.{w.af}_othaf.{w.othaf}_cov.{w.cov}_othcov.{w.othcov}_ncells.{w.ncells}_othncells.{w.othncells}_mean.{w.mean}/output.ipynb"
    return



def get_repr_cl_note(wildcards):
    w = wildcards
    if len(condition) == 1:
        return join(ROOT_DIR,  "workflow/notebooks/individual_clones/oneCondition_get_representative_clones.ipynb"),
    return join(ROOT_DIR,  "workflow/notebooks/individual_clones/get_representative_clones.ipynb"),


def get_rank_cl_sizes_condition(wildcards):
    if len(condition) == 1:
        return condition[0]
    return "None"


rule get_repr_clones_cl:
    input:
        hyper = get_hypergeo,
        cells="{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/cells_meta.tsv"
    output:
        note = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/clones_ranked/output_rank_ncells.ipynb",
        clone_order = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/clones_ranked/representative_cloneID.txt"
    params:
        indir = lambda wildcards, input: dirname(input.hyper),
        outdir = lambda wildcards, output: dirname(output.note),
        note = get_repr_cl_note, #note = join(ROOT_DIR,  "workflow/notebooks/individual_clones/rank_clone_lineage_sizes.ipynb"),
        is_mt = False,
        condition = get_rank_cl_sizes_condition

    shell: "papermill -p indir {params.indir} -p outdir {params.outdir} -p cells_meta_f {input.cells} -p is_mt {params.is_mt} -p condition {params.condition} {params.note} {output.note}"


rule get_repr_clones_mt:
    input:
        hyper = get_hypergeo,
        cells = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{clone_shift_method}/cells_meta.tsv"
    output:
        note = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{clone_shift_method}/clones_ranked/output_rank_ncells.ipynb",
        clone_order = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{clone_shift_method}/clones_ranked/representative_cloneID.txt"
    params:
        indir = lambda wildcards, input: dirname(input.hyper),
        outdir = lambda wildcards, output: dirname(output.note),
        note = get_repr_cl_note,
        is_mt = lambda wildcards: True if wildcards.clone_shift_method=="mt_as_clones" else False,
        condition = get_rank_cl_sizes_condition
        #p_thresh = 0.1
    shell: "papermill -p indir {params.indir} -p outdir {params.outdir} -p cells_meta_f {input.cells} -p is_mt {params.is_mt} -p condition {params.condition} {params.note} {output.note}"


rule preproc_complete:
    input:
        expand("{{outdir}}/{{clone_subset}}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_mt_as_clones/cells_meta.tsv",
                d=np.arange(config["N_DONORS"]),
                variants=[x for x in params_clones["variants"] if x != "simple"],
                af= best_p["af"], othaf= best_p["oth_af"],
                cov= best_p["cov"], othcov= best_p["oth_cov"],
                ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
                mean= best_p["mean_pos_cov"]),
        expand("{{outdir}}/{{clone_subset}}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/{outputs}",
                  outputs=["clones_ranked/output_rank_ncells.ipynb"], #"cloneIDs.txt","clones_ranked/output.ipynb",
                  d=np.arange(config["N_DONORS"]), variants=[x for x in params_clones["variants"] if x != "simple"],
                  clone_shift_method=knn_clone_shift,
                  kparam=params_clones["knn"]["params"]["resolution"]),
        expand("{{outdir}}/{{clone_subset}}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{clone_shift_method}/{outputs}",
                  outputs=["clones_ranked/output_rank_ncells.ipynb"], #"cloneIDs.txt","clones_ranked/output.ipynb"
                  d=np.arange(config["N_DONORS"]), clone_shift_method=mt_clone_shift,
                  variants=[x for x in params_clones["variants"] if x != "simple"],
                  af= best_p["af"], othaf= best_p["oth_af"],
                  cov= best_p["cov"], othcov= best_p["oth_cov"],
                  ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
                  mean= best_p["mean_pos_cov"],),
    output: "{outdir}/{clone_subset}/.preproc.txt"
    shell: "touch {output}"


checkpoint get_single_cloneID:
    input:
        clone_f = "{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/cloneIDs.txt",
    output:
        clone_dir = directory("{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/cloneIDs")
    run:
        if not os.path.exists(output.clone_dir):
            os.mkdir(output.clone_dir)
        with open(input.clone_f, "r") as f:
            lines = f.readlines()
        lines = [x.strip() for x in lines]
        ic('lines', lines)
        for c_id in lines:
            ic('c_id', c_id)
            with open(join(output.clone_dir, f"{c_id}.txt")) as f:
                f.write(str(c_id))


# input function for rule aggregate, return paths to all files produced by the checkpoint 'somestep'
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.get_single_cloneID.get(**wildcards).output[0]
    return expand("{{outdir}}/{{clone_subset}}/donor{{d}}/cloneMethod_{{method}}/clonalShift_method_{{clone_shift_method}}/single/{cloneID}/hypergeo_pval.png",
           cloneID=glob_wildcards(os.path.join(checkpoint_output, "{cloneID}.txt")).i)


# an aggregation over all produced clusters
rule aggregate:
    input:
        aggregate_input
    output:
        "{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/aggregated.txt"
    shell:
        "touch {output}"



""" Get clone_overlay of UMAP. 
"""
def get_se(wildcards):
    w = wildcards
    clshift = w.clone_shift_method
    if  clshift == "dendro_bc" or clshift == "clones":
        return f"clones/variants_{w.variants}/knn/kparam_{w.kparam}/gff_{config['gff']}/annotation_clones/SE.rds"
        #return f"{w.outdir}/anno_multiplex/gff_{config['gff']}/se.rds"
    elif clshift == "mt_as_clones_dendro" or clshift == "mt_as_clones":
        return f"clones/variants_{w.variants}/knn/kparam_30/gff_{config['gff']}/annotation_clones/SE.rds"
        #return f"{w.outdir}/anno_multiplex/gff_{config['gff']}/se.rds"

def get_mt_f(wildcards):
    w = wildcards
    outdir = join(w.outdir, f"multiplex/clones_{w.variants}/donor{w.d}/af.tsv")
    return outdir


def get_top_clone_hypergeoSig_script(wildcards):
    return join(ROOT_DIR, "workflow/notebooks/individual_clones/oneCondition_top_individual_clone_lineage_Sig_hypergeo.ipynb")


def get_param_top_clone_hypergeoSig_condition():
    if len(condition) == 1:
        return condition[0]
    else:
        return "None"


rule top_clone_hypergeoSig_cl:
    input:
        cloneIDs_f = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/clones_ranked/representative_cloneID.txt",
        indir = get_hypergeo,
    output:
        note = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/top/top_hypergeo_sig.ipynb",
    params:
        notebook = get_top_clone_hypergeoSig_script, #join(ROOT_DIR, "workflow/notebooks/individual_clones/top_individual_clone_lineage_Sig_hypergeo.ipynb"),
        indir = lambda wildcards, input: dirname(input.indir),
        outdir = lambda wildcards, output: dirname(output[0]),
        condition=get_param_top_clone_hypergeoSig_condition(),
        ntop_clones = 10
        #clone_id = lambda wildcards: wildcards.cloneID
    shell:
        "papermill -p indir {params.indir} -p outdir {params.outdir} -p cloneID_f {input.cloneIDs_f} -p condition {params.condition} -p ntop_clones {params.ntop_clones} {params.notebook} {output.note}"


rule top_clone_hypergeoSig_mt:
    input:
        cloneIDs_f = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{clone_shift_method}/clones_ranked/representative_cloneID.txt",
        indir = get_hypergeo,
    output:
        note = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{clone_shift_method}/top/top_hypergeo_sig.ipynb",
    params:
        notebook = get_top_clone_hypergeoSig_script,
        indir = lambda wildcards, input: dirname(input.indir),
        outdir = lambda wildcards, output: dirname(output[0]),
        p_thresh = 0.1,
        ntop_clones = 10,
        condition=get_param_top_clone_hypergeoSig_condition()
        #clone_id = lambda wildcards: wildcards.cloneID
    shell:
        "papermill -p indir {params.indir} -p outdir {params.outdir} -p cloneID_f {input.cloneIDs_f} -p ntop_clones {params.ntop_clones} -p condition {params.condition} -p p_thresh {params.p_thresh} {params.notebook} {output.note}"


rule top_clone_lineage_count:
    input:
        cells = "{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/cells_meta.tsv",
        clone_order_f = "{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/clones_ranked/representative_cloneID.txt",
        #se =  expand("{{outdir}}/anno_multiplex/gff_{gff}/SE.rds", gff=config["gff"])
    output:
        note = "{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/top/top_clone_lineage_count.ipynb",
        #clust_cond_f = "{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/top/top_cluster_condition_ncells.pdf"
    params:
        script = join(ROOT_DIR, "workflow/notebooks/individual_clones/top_individual_clone_lineage_distribution.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
        ntop_clones = 10
    shell: "papermill -p cells_meta_f {input.cells} -p outdir {params.outdir} -p clone_order_f {input.clone_order_f} -p ntop_clones {params.ntop_clones} {params.script} {output.note}"


def get_clone_mt_variants_script(wildcards):
    w = wildcards
    if w.clone_order == "top":
        return join(ROOT_DIR, "workflow/notebooks/individual_clones/top_individual_clone_mt_variants.ipynb"),
    elif w.clone_order == "all":
        return join(ROOT_DIR, "workflow/notebooks/individual_clones/all_individual_clone_mt_variants.ipynb"),
    else:
        raise ValueError


rule top_clone_mt_variants_cl:
    input:
        cells = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{cloneShift_method}/cells_meta.tsv",
        clone_order_f = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{cloneShift_method}/clones_ranked/representative_cloneID.txt",
        mt_f = get_mt_f
    output:
        note = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{cloneShift_method}/top/{clone_order}_clone_mt_variants.ipynb",
    params:
        indir = lambda wildcards, input: dirname(input.mt_f),
        outdir = lambda wildcards, output: dirname(output.note),
        script = get_clone_mt_variants_script, #join(ROOT_DIR, "workflow/notebooks/individual_clones/top_individual_clone_mt_variants.ipynb"),
        af_thresh=0.001,
        cov_thresh=2,
        ntop_clones=10,
        ntop_vars=10
    shell: "papermill -p cells_meta_f {input.cells} -p indir {params.indir} -p outdir {params.outdir} -p clone_order_f {input.clone_order_f} -p ntop_clones {params.ntop_clones} -p af_thresh {params.af_thresh} -p cov_thresh {params.cov_thresh} -p ntop_vars {params.ntop_vars} {params.script} {output.note}"


rule top_clone_mt_variants_mt:
    input:
        cells = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{cloneShift_method}/cells_meta.tsv",
        clone_order_f = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{cloneShift_method}/clones_ranked/representative_cloneID.txt",
        mt_f = get_mt_f
    output:
        note = "{outdir}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{cloneShift_method}/top/{clone_order}_clone_mt_variants.ipynb",
    params:
        indir = lambda wildcards, input: dirname(input.mt_f),
        outdir = lambda wildcards, output: dirname(output.note),
        script = get_clone_mt_variants_script, #join(ROOT_DIR, "workflow/notebooks/individual_clones/top_individual_clone_mt_variants.ipynb"),
        af_thresh=0.001,
        cov_thresh=2,
        ntop_clones=20,
        ntop_vars=10
    shell: "papermill -p cells_meta_f {input.cells} -p indir {params.indir} -p outdir {params.outdir} -p clone_order_f {input.clone_order_f} -p ntop_clones {params.ntop_clones} -p af_thresh {params.af_thresh} -p cov_thresh {params.cov_thresh} -p ntop_vars {params.ntop_vars} {params.script} {output.note}"



def get_top_clone_umap_overlay_script(wildcards):
    if wildcards.cloneShift_method == "mt_as_clones":
        return join(ROOT_DIR, "workflow/notebooks/individual_clones/mt_top_individual_clones_umap_overlap.ipynb")
    else:
        return join(ROOT_DIR, "workflow/notebooks/individual_clones/top_individual_clones_umap_overlap.ipynb")


rule top_clone_umap_overlay:
    input:
        cells = "{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/cells_meta.tsv",
        clone_order_f = "{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/clones_ranked/representative_cloneID.txt",
        se =  expand("{{outdir}}/anno_multiplex/gff_{gff}/SE.rds", gff=config["gff"])
    output:
        note = "{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/top/top_umap_overlay.ipynb"
    params:
        script = get_top_clone_umap_overlay_script, #join(ROOT_DIR, "workflow/notebooks/individual_clones/top_individual_clones_umap_overlap.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
    shell: "papermill -p se_f {input.se} -p cells_meta_f {input.cells} -p outdir {params.outdir} -p clone_order_f {input.clone_order_f} {params.script} {output.note}"


rule top_merge:
    input:
        ins = multiext("{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/top/","top_umap_overlay.ipynb", "top_hypergeo_sig.ipynb", "top_clone_lineage_count.ipynb", "top_clone_mt_variants.ipynb"),
    output:
        out_f = multiext("{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/top/","clone_shift_combine.pdf"),
        note = "{outdir}/{clone_subset}/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/top/merge_panels.ipynb"
    params:
        outdir = lambda wildcards, output:  dirname(output.note),
        script = join(ROOT_DIR, "workflow/notebooks/individual_clones/concat_images_top_clones.ipynb"),
    shell: "papermill -p outdir {params.outdir} {params.script} {output.note}"


rule top_clone_complete:
    input:
        expand("{{outdir}}/{{clone_subset}}/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/top/{f}",
              f=["top_umap_overlay.ipynb", "top_hypergeo_sig.ipynb", "top_clone_lineage_count.ipynb", "top_clone_mt_variants.ipynb", "all_clone_mt_variants.ipynb", "clone_shift_combine.pdf"],
              d=np.arange(config["N_DONORS"]), variants=[x for x in params_clones["variants"] if x != "simple"],
              clone_shift_method = knn_clone_shift,
              kparam=params_clones["knn"]["params"]["resolution"]),

        # expand("{{outdir}}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{clone_shift_method}/top/{f}",
        #           f=["top_umap_overlay.ipynb", "top_hypergeo_sig.ipynb", "top_clone_lineage_count.ipynb", "top_clone_mt_variants.ipynb", "all_clone_mt_variants.ipynb", "clone_shift_combine.pdf"],
        #           d=np.arange(config["N_DONORS"]), clone_shift_method=mt_clone_shift,
        #           variants=[x for x in params_clones["variants"] if x != "simple"],
        #           af= best_p["af"], othaf= best_p["oth_af"],
        #           cov= best_p["cov"], othcov= best_p["oth_cov"],
        #           ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
        #           mean= best_p["mean_pos_cov"]),
        # expand("{{outdir}}/{clone_subset}/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_mt_as_clones/top/{f}",
        #           f=["top_clone_lineage_count.ipynb", "top_hypergeo_sig.ipynb", "top_umap_overlay.ipynb"],
        #           d=np.arange(config["N_DONORS"]), clone_shift_method=mt_clone_shift,
        #           variants=[x for x in params_clones["variants"] if x != "simple"],
        #           af= best_p["af"], othaf= best_p["oth_af"],
        #           cov= best_p["cov"], othcov= best_p["oth_cov"],
        #           ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
        #           mean= best_p["mean_pos_cov"]),
    output: "{outdir}/{clone_subset}/.top_clones.txt"
    shell: "touch {output}"