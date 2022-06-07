import pandas as pd
from src.config import ROOT_DIR
from src.utils.parse_config import read_config_file
import os
import numpy as np
from os.path import join, dirname
from snakemake.utils import min_version
from icecream import ic
min_version("6.0")



knn_clone_shift= ["clones", "dendro_bc"]
#mt_method = "mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}"
mt_clone_shift = ["mt_as_clones_dendro"]#, "mt_as_clones"]
best_p = config["mt_as_clones"]["best_params"]
params_clones = config["clones"]


dendro_d = config['clones']['params']['dendro_thresh']

def load_cells_meta(cells_meta_f):
    print(cells_meta_f)
    cells_meta = pd.read_csv(cells_meta_f, sep="\t", index_col=0)
    cells_meta = cells_meta.loc[~(cells_meta["name"]=="None")]
    if not "cluster_labels" in cells_meta.columns.values:
        cells_meta["clusterID"] = cells_meta["seurat_clusters"]
    else:
        cells_meta["clusterID"] = cells_meta["cluster_labels"]
    return cells_meta




""" Output line-separated cell IDs
Cell ID, Clone, Lineage, Condition
"""
rule get_clone_cells:
    input:
        se_cells_meta_f = "{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/gff_A2_black/annotation_clones/se_cells_meta_labels.tsv"
    output: "{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_clones/cells_meta.tsv"
    run:
        # cells_meta = pd.read_csv(input.se_cells_meta_f, sep="\t", index_col=0)
        # cells_meta = cells_meta.loc[~(cells_meta["name"]=="None")]
        # if not "cluster_labels" in cells_meta.columns.values:
        #     cells_meta["clusterID"] = cells_meta["seurat_clusters"]
        # else:
        #     cells_meta["clusterID"] = cells_meta["cluster_labels"]
        cells_meta = load_cells_meta(input.se_cells_meta_f)
        cells_meta = cells_meta.rename({"name":"cloneID"}, axis=1)
        cells_meta = cells_meta.loc[cells_meta["donor"] == wildcards.d]
        cells_meta[["cloneID", "clusterID","condition", "donor"]].to_csv(output[0], sep="\t")


rule get_clone_dendro_cells:
    input:
        barcodes_dir = expand("{{outdir}}/clones/variants_{{variants}}/knn/kparam_{{kparam}}/barcodes/btwnClones_dendro_dt_{dt}/donor{{d}}.clones_dendro.csv", dt=dendro_d),
        se_cells_meta_f = "{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/gff_A2_black/annotation_clones/se_cells_meta_labels.tsv"
    output:
        cells = "{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_dendro_bc/cells_meta.tsv"
    #output: "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/{clone}/cells.txt"
    run:
        clone_col = "den_clust"
        cells_meta = load_cells_meta(input.se_cells_meta_f)
        #barcodes_in = pd.read_csv(join(input.barcodes_dir[0],f"donor{wildcards.d}.clones_dendro.csv"), index_col=0)
        barcodes_in = pd.read_csv(list(input.barcodes_dir)[0], index_col=0)
        print('d', wildcards.d)
        #print('barcodes_in')
        print('a',barcodes_in.head())
        barcodes_in[clone_col] = str(wildcards.d) + "_" + barcodes_in[clone_col]
        print('b',barcodes_in.head())
        cells_meta = cells_meta[cells_meta["name"].isin(barcodes_in.index)]
        cells_meta[clone_col] = cells_meta.apply(lambda x: barcodes_in.loc[x["name"], clone_col] , axis=1)
        print('c', cells_meta.head())
        cells_meta = cells_meta.loc[cells_meta["donor"].astype(str) == str(wildcards.d)]
        cells_meta = cells_meta.rename({clone_col: "cloneID"}, axis=1)
        cells_meta[["cloneID", "clusterID","condition", "donor"]].to_csv(output.cells, sep="\t")


rule get_cells_mt_as_clones:
    input:
        clones = "{outdir}/mt_as_clones/variants_{variants}/bestparams/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/best_params_save.ipynb",
        se_cells_meta_f = expand("{{outdir}}/anno_multiplex/gff_{gff}/se_cells_meta_labels.tsv", gff=config["gff"])
    output:
        cells="{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_mt_as_clones/cells_meta.tsv"
        #"{outdir}/single_clones/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/{clone}/cells.txt"
    params:
        indir = lambda wildcards, input: dirname(input.clones)
    run:
        indir = params.indir
        best_params = pd.read_csv(join(indir,  "best_params.csv"))
        bin_d = pd.read_csv(join(indir,  f"donor_{wildcards.d}_binary.csv"), index_col=0)
        bin_d.columns = [f"cloneID_{x}" for x in bin_d.columns]
        labels_df = pd.read_csv(list(input.se_cells_meta_f)[0],sep="\t").set_index("ID")
        print(labels_df.shape)
        labels_df = labels_df[~(labels_df["donor"]=='None')]
        print(labels_df.shape)
        cells_meta = pd.merge(labels_df[["donor", "cluster_labels", "condition"]], bin_d, left_index=True, right_index=True, how="inner")
        cells_meta.index = [f'{x.split("_")[1]}_{x.split("_")[0]}'  for x in cells_meta.index]
        cells_meta = cells_meta.rename({"cluster_labels":"clusterID"},axis=1)
        cells_meta.to_csv(output.cells, sep="\t")

# convert back to raw cell count from log2
#lin_mt_ncells["count"] = (np.ceil((2**lin_mt_ncells["count"])-1)).astype(int)

rule get_cells_mt_as_clones_dendro:
    input:
        clones = "{outdir}/mt_as_clones/variants_{variants}/dendro/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/dendro_mt_clones.ipynb",
        se_cells_meta_f = expand("{{outdir}}/anno_multiplex/gff_{gff}/se_cells_meta_labels.tsv", gff=config["gff"])
    output:
        cells="{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_mt_as_clones_dendro/cells_meta.tsv"
    params:
        indir = lambda wildcards, input: dirname(input.clones)
    run:
        indir = params.indir
        clone_col="den_clust"
        atac_col = "cluster_labels"
        labels_df = pd.read_csv(list(input.se_cells_meta_f)[0],sep="\t").reset_index().set_index("ID")
        print(labels_df.shape)
        labels_df = labels_df[~(labels_df["donor"]=='None')]
        print(labels_df.shape)
        labels_df.head()
        den_d = pd.read_csv(join(indir, f"don_{wildcards.d}_mt_dendro_clust.csv"), index_col=0)
        labels = labels_df[labels_df["donor"] == wildcards.d]
        cells_meta = pd.merge(left=den_d, right=labels[[atac_col]], left_index=True, right_index=True, how="inner")
        cells_meta.index = [f'{x.split("_")[1]}_{x.split("_")[0]}'  for x in cells_meta.index]
        cells_meta = cells_meta.rename({"Donor":"donor", "cluster_labels":"clusterID", clone_col:"cloneID"},axis=1)
        cells_meta.to_csv(output.cells, sep="\t")


rule get_cloneIDs:
    input:
        cells = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/cells_meta.tsv"
    output:
        clone_f = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/cloneIDs.txt",
    run:
        df = pd.read_csv(input.cells, sep="\t")
        if wildcards.clone_shift_method == "mt_as_clones":
            cloneIDs = [x for x in df.columns if "cloneID" in x ]
        else:
            cloneIDs = list(set(df["cloneID"]))
        print('cloneIDs', cloneIDs)

        cloneIDs_str = "\n".join(cloneIDs)
        with open(output.clone_f, "w") as f:
            f.write(cloneIDs_str)

        # for c_id, val in enumerate(cloneIDs):
        #     with open(join(outdir, f"{c_id}.txt")):
        #         f.write(str(c_id))


# input function for rule aggregate, return paths to all files produced by the checkpoint 'somestep'
# def aggregate_input(wildcards):
#     checkpoint_output = checkpoints.get_cloneIDs.get(**wildcards).output[0]
#     with open(checkpoint_output, 'r') as f:
#         lines = f.readlines()
#     lines = [x.strip() for x in lines]
#     print('lines first and list', lines[0], lines[-1])
#     return
#     return expand("my_directory/{i}.txt",
#                   i=glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)


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
#clonal_shifts/donors/donor0/variants_init/knn/kparam_30/clones/"


rule get_rank_mt:
    input:
        hyper = get_hypergeo,
        #indir = "/data/Mito_Trace/output/pipeline/v02/CHIP_b1/MTBlacklist_A2/data/merged/MT/cellr_True/numread_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/mgatk/vireoIn/clonal_shifts/variants_init/donors/donor0/clones/knn_kparam_30/"
    output:
        #note = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/clones_ranked/output.ipynb",
        note = "{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{clone_shift_method}/clones_ranked/output.ipynb",
    params:
        indir = lambda wildcards, input: dirname(input.hyper),
        outdir = lambda wildcards, output: dirname(output.note),
        note = join(ROOT_DIR,  "workflow/notebooks/individual_clones/rank_clone_lineage_hypergeo.ipynb"),
        #p_thresh = 0.1
    shell: "papermill -p indir {params.indir} -p outdir {params.outdir} {params.note} {output.note}"



rule get_rank_cl:
    input:
        hyper = get_hypergeo
    output:
        note = "{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/clones_ranked/output.ipynb",
    params:
        indir = lambda wildcards, input: dirname(input.hyper),
        outdir = lambda wildcards, output: dirname(output.note),
        note = join(ROOT_DIR,  "workflow/notebooks/individual_clones/rank_clone_lineage_hypergeo.ipynb"),
        #p_thresh = 0.1
    shell: "papermill -p indir {params.indir} -p outdir {params.outdir} {params.note} {output.note}"




#/data/Mito_Trace/output/pipeline/v02/CHIP_b1/MTBlacklist_A2/data/merged/MT/cellr_True/numread_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/mgatk/vireoIn/single_clones/donor0/cloneMethod_variants_init_knn_resolution_3/clonalShift_method_dendro_bc/cells_meta.tsv
rule preperoc_complete:
    input:
        expand("{{outdir}}/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/{outputs}",
                  outputs=["cloneIDs.txt","clones_ranked/output.ipynb"],
                  d=np.arange(config["N_DONORS"]), variants=[x for x in params_clones["variants"] if x != "simple"],
                  clone_shift_method=knn_clone_shift,
                  kparam=params_clones["knn"]["params"]["resolution"]),
        expand("{{outdir}}/single_clones/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{clone_shift_method}/{outputs}",
                  outputs=["cloneIDs.txt","clones_ranked/output.ipynb"],
                  d=np.arange(config["N_DONORS"]), clone_shift_method=mt_clone_shift,
                  variants=[x for x in params_clones["variants"] if x != "simple"],
                  af= best_p["af"], othaf= best_p["oth_af"],
                  cov= best_p["cov"], othcov= best_p["oth_cov"],
                  ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
                  mean= best_p["mean_pos_cov"],),
    output: "{outdir}/single_clones/.preproc.txt"
    shell: "touch {output}"


# the checkpoint that shall trigger re-evaluation of the DAG
# an number of file is created in a defined directory
# checkpoint somestep:
#     output:
#         directory("my_directory/")
#     shell:
#         "mkdir my_directory/;"
#         "for i in 1 2 3; do touch $i.txt; done"


checkpoint get_single_cloneID:
    input:
        clone_f = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/cloneIDs.txt",
    output:
        clone_dir = directory("{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/cloneIDs")
    run:
        if not os.path.exists(output.clone_dir):
            os.mkdir(output.clone_dir)
        with open(input.clone_f, "r") as f:
            lines = f.readlines()
        lines = [x.strip() for x in lines]
        print('lines', lines)
        for c_id in lines:
            print('c_id', c_id)
            with open(join(output.clone_dir, f"{c_id}.txt")) as f:
                f.write(str(c_id))


def get_rank():
    return

#cloneID = ""


# rule ind_clone_hypergeo:
#     input:
#         cells = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/cells_meta.tsv",
#         indir = get_hypergeo,
#         cloneID = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/cloneIDs/{cloneID}.txt"
#     output:
#         note = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/single/{cloneID}/hypergeo_pval.ipynb",
#         fig = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/single/{cloneID}/hypergeo_pval.png"
#         # multiext("{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/{cloneID}/",
#         #               "hypergeo_pval.png", "clone_shift.png")
#     params:
#         notebook = join(ROOT_DIR, "workflow/notebooks/individual_clones/individual_clone_lineage_hypergeo.ipynb"),
#         indir = lambda wildcards, input: dirname(input.indir),
#         outdir = lambda wildcards, output: dirname(output[0]),
#         clone_id = lambda wildcards: wildcards.cloneID
#     shell:
#         "papermill -p indir {input.indir} -p outdir {params.outdir} -p clone_id {params.clone_id} {params.notebook} {output.note}"


# input function for rule aggregate, return paths to all files produced by the checkpoint 'somestep'
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.get_single_cloneID.get(**wildcards).output[0]
    return expand("{{outdir}}/single_clones/donor{{d}}/cloneMethod_{{method}}/clonalShift_method_{{clone_shift_method}}/single/{cloneID}/hypergeo_pval.png",
           cloneID=glob_wildcards(os.path.join(checkpoint_output, "{cloneID}.txt")).i)


# #rule run_clone_umap:
# rule aggregate:
#     input:
#         dynamic("{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/{cloneID}/hypergeo_pval.png")
#     output:
#         "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/aggegated.txt"
#     shell:
#         "touch {output}"


# an aggregation over all produced clusters
rule aggregate:
    input:
        aggregate_input
    output:
        "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/aggregated.txt"
    shell:
        "touch {output}"



""" Get clone_overlay of UMAP. 
"""
# rule ind_clone_umap_overlay:
#     input:  "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/{clone}/cells_meta.tsv"
#     output:
#         multiext("{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/{clone}/{clone}.", "umap_clust_color.png", "umap_cond_color.png")
#     params:
#         script = join(ROOT_DIR, "workflow/notebooks/individual_clones/individual_clones_umap_overlap.ipynb")

def get_se(wildcards):
    w = wildcards
    clshift = w.clone_shift_method
    if  clshift == "dendro_bc" or clshift == "clones":
        return f"clones/variants_{w.variants}/knn/kparam_{w.kparam}/gff_{config['gff']}/annotation_clones/SE.rds"
        #return f"{w.outdir}/anno_multiplex/gff_{config['gff']}/se.rds"
    elif clshift == "mt_as_clones_dendro" or clshift == "mt_as_clones":
        return f"clones/variants_{w.variants}/knn/kparam_30/gff_{config['gff']}/annotation_clones/SE.rds"
        #return f"{w.outdir}/anno_multiplex/gff_{config['gff']}/se.rds"


rule ind_clone_hypergeo_cl:
    input:
        cloneIDs_f = "{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/cloneIDs.txt",
        indir = get_hypergeo,
        #cloneID = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/cloneIDs/{cloneID}.txt"
    output:
        note = "{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/single/hypergeo_pval.ipynb",
        #fig = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/single/{cloneID}/hypergeo_pval.png"
        # multiext("{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/{cloneID}/",
        #               "hypergeo_pval.png", "clone_shift.png")
    params:
        notebook = join(ROOT_DIR, "workflow/notebooks/individual_clones/individual_clone_lineage_hypergeo.ipynb"),
        indir = lambda wildcards, input: dirname(input.indir),
        outdir = lambda wildcards, output: dirname(output[0]),
        #clone_id = lambda wildcards: wildcards.cloneID
    shell:
        "papermill -p indir {params.indir} -p outdir {params.outdir} -p cloneID_f {input.cloneIDs_f} {params.notebook} {output.note}"


rule ind_clone_hypergeo_mt:
    input:
        cloneIDs_f = "{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{clone_shift_method}/cloneIDs.txt",
        indir = get_hypergeo,
        #cloneID = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/cloneIDs/{cloneID}.txt"
    output:
        note = "{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{clone_shift_method}/single/hypergeo_pval.ipynb",
        #fig = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{clone_shift_method}/single/{cloneID}/hypergeo_pval.png"
        # multiext("{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/{cloneID}/",
        #               "hypergeo_pval.png", "clone_shift.png")
    params:
        notebook = join(ROOT_DIR, "workflow/notebooks/individual_clones/individual_clone_lineage_hypergeo.ipynb"),
        indir = lambda wildcards, input: dirname(input.indir),
        outdir = lambda wildcards, output: dirname(output[0]),
        #clone_id = lambda wildcards: wildcards.cloneID
    shell:
        "papermill -p indir {params.indir} -p outdir {params.outdir} -p cloneID_f {input.cloneIDs_f} {params.notebook} {output.note}"


rule ind_clone_umap_overlay:
    input:
        cells = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/cells_meta.tsv",
        se =  expand("{{outdir}}/anno_multiplex/gff_{gff}/SE.rds", gff=config["gff"])
    output:
        note = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/single/umap_overlay.ipynb"
    params:
        script = join(ROOT_DIR, "workflow/notebooks/individual_clones/individual_clones_umap_overlap.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note)
    shell: "papermill -p se_f {input.se} -p cells_meta_f {input.cells} -p outdir {params.outdir}  {params.script} {output.note}"


rule ind_clone_lineage_count:
    input:
        cells = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/cells_meta.tsv",
        #se =  expand("{{outdir}}/anno_multiplex/gff_{gff}/SE.rds", gff=config["gff"])
    output:
        note = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/single/clone_lineage_count.ipynb"
    params:
        script = join(ROOT_DIR, "workflow/notebooks/individual_clones/individual_clone_lineage_distribution.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note)
    shell: "papermill -p cells_meta_f {input.cells} -p outdir {params.outdir} {params.script} {output.note}"


def get_mt_f(wildcards):
    w = wildcards
    outdir = join(w.outdir, f"multiplex/clones_{w.variants}/donor{w.d}/af.tsv")
    return outdir


# rule ind_clone_coverage_af:
#     input:
#         cells = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/cells_meta.tsv",
#         mt_indir = get_mt_indir
#     output:
#         note = "{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/single/mt_coverage.ipynb"
#     params:
#         script = join(ROOT_DIR, "workflow/notebooks/individual_clones/individual_clone_mt_coverage.ipynb"),
#         outdir = lambda wildcards, output: dirname(output.note)
#     shell: "papermill -p indir {input.mt_indir} -p cells_meta_f {input.cells} -p outdir {params.outdir} -d DONOR {wildcards.d {params.script} {output.note}"

#multiplex/clones_
rule ind_clone_coverage_af_cl:
    input:
        cells = "{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{cloneShift_method}/cells_meta.tsv",
        mt_f = get_mt_f
    output:
        note = "{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{cloneShift_method}/single/mt_coverage.ipynb"
    params:
        mt_indir = lambda wildcards, input: dirname(input.mt_f),
        script = join(ROOT_DIR, "workflow/notebooks/individual_clones/individual_clone_mt_coverage.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note)
    shell: "papermill -p indir {params.mt_indir} -p cells_meta_f {input.cells} -p outdir {params.outdir} {params.script} {output.note}"


rule ind_clone_coverage_af_mt:
    input:
        cells = "{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{cloneShift_method}/cells_meta.tsv",
        mt_f = get_mt_f
    output:
        note = "{outdir}/single_clones/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{cloneShift_method}/single/mt_coverage.ipynb"
    params:
        mt_indir = lambda wildcards, input: dirname(input.mt_f),
        script = join(ROOT_DIR, "workflow/notebooks/individual_clones/individual_clone_mt_coverage.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note)
    shell: "papermill -p indir {params.mt_indir} -p cells_meta_f {input.cells} -p outdir {params.outdir} {params.script} {output.note}"



# """get distinct variants"""
# rule get_variants:
#     input: ""
#     output: ""


rule individual_clone_complete:
    input:
        #expand("{outdir}/single_clones/donor{d}/cloneMethod_{method}/clonalShift_method_{cloneShift_method}/aggegated.txt"),
        # expand("{{outdir}}/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/aggregated.txt",
        #       d=np.arange(config["N_DONORS"]), variants=[x for x in params_clones["variants"] if x != "simple"],
        #       clone_shift_method = knn_clone_shift,
        #       kparam=params_clones["knn"]["params"]["resolution"]),
        # expand("{{outdir}}/single_clones/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{clone_shift_method}/aggregated.txt",
        #           d=np.arange(config["N_DONORS"]), clone_shift_method=mt_clone_shift,
        #           variants=[x for x in params_clones["variants"] if x != "simple"],
        #           af= best_p["af"], othaf= best_p["oth_af"],
        #           cov= best_p["cov"], othcov= best_p["oth_cov"],
        #           ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
        #           mean= best_p["mean_pos_cov"]),
        expand("{{outdir}}/single_clones/donor{d}/cloneMethod_variants_{variants}_knn_resolution_{kparam}/clonalShift_method_{clone_shift_method}/single/{f}",
              f=["umap_overlay.ipynb", "hypergeo_pval.ipynb", "clone_lineage_count.ipynb", "mt_coverage.ipynb"],
              d=np.arange(config["N_DONORS"]), variants=[x for x in params_clones["variants"] if x != "simple"],
              clone_shift_method = knn_clone_shift,
              kparam=params_clones["knn"]["params"]["resolution"]),
        expand("{{outdir}}/single_clones/donor{d}/cloneMethod_variants_{variants}_mt_bestparams_af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/clonalShift_method_{clone_shift_method}/single/{f}",
                  f=["umap_overlay.ipynb", "hypergeo_pval.ipynb", "clone_lineage_count.ipynb", "mt_coverage.ipynb"],
                  d=np.arange(config["N_DONORS"]), clone_shift_method=mt_clone_shift,
                  variants=[x for x in params_clones["variants"] if x != "simple"],
                  af= best_p["af"], othaf= best_p["oth_af"],
                  cov= best_p["cov"], othcov= best_p["oth_cov"],
                  ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
                  mean= best_p["mean_pos_cov"]),
    output: "{outdir}/single_clones/.aggregate.txt"
    shell: "touch {output}"