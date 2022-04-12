from os.path import join, dirname
from src.utils.data_io import sparse_to_cellranger_fragments
from src.config import ROOT_DIR

def get_anno_integrate(wildcards):
    print(join(config["anno_res"], "mergedSamples","allSamples.integrated.rds"))
    return join(config["anno_res"], "mergedSamples","allSamples.integrated.rds")

params = config["annotation_clones"]["params"]
params_cl_embed = config["clone_clust_embed"]["params"]

rule addClones:
    input:
        noc = get_anno_integrate, #"{outdir}/annotation_clones/mergedSamples/allSamples.integrated.rds", "{anno_indir}/mergedSamples/allSamples.integrated.rds", #
        clones = "{outdir}/cells_meta.tsv",#"{outdir}/pipeline/published/{prefix}/data/clones/clones.txt",
    output:
        se_f = "{outdir}/annotation_clones/SE.rds",
        note = "{outdir}/annotation_clones/addClones.ipynb"
    params:
        outdir = lambda wildcards, output: dirname(output[0]),
        rscript= join(ROOT_DIR, "workflow/notebooks/lineage_clones/addClones.ipynb"), # The script defaults to the granja data
    shell: "papermill -p cells_meta_f {input.clones} -p se_f {input.noc} -p outdir {params.outdir} {params.rscript} {output.note}"


rule se_meta:
    """Prepare clone-by-cluster counts for umap and hypergeometric test"""
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
    output:
        se_meta = "{outdir}/annotation_clones/se_cells_meta.tsv",
        note = "{outdir}/annotation_clones/se_cells_meta.ipynb",
    params:
        rscript = join(ROOT_DIR, "workflow/notebooks/lineage_clones/clusters_cells_meta.ipynb"),
    shell: "papermill -p se_f {input.se_f} {params.rscript} {output.note}"

def get_cluster_labels():
    return config.get("umap_clusters_f", "FALSE")


rule add_cluster_labels:
    """Prepare clone-by-cluster counts for umap and hypergeometric test"""
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
    output:
        se_meta = "{outdir}/annotation_clones/se_cells_meta_labels.tsv",
        note = "{outdir}/annotation_clones/add_cluster_labels.ipynb",
    params:
        labels = get_cluster_labels(),
        rscript = join(ROOT_DIR, "workflow/notebooks/lineage_clones/add_cluster_labels.ipynb"),
    shell: "papermill -p se_f {input.se_f} -p cluster_labels_f {params.labels} {params.rscript} {output.note}"


def get_lineage_clone_counts_script(wildcards):
    if wildcards.donType == "combinedDonors":
        return join(ROOT_DIR, "workflow/notebooks/lineage_clones/clone_cluster_input.ipynb")
    return join(ROOT_DIR, "workflow/notebooks/lineage_clones/clone_cluster_input_donors.ipynb")

# Barplots of lineage and clones
rule lineage_clone_counts:
    """ ISSUE with the output the size is not accurate"""
    input:
        se_meta = "{outdir}/annotation_clones/se_cells_meta_labels.tsv",
    output:
        note = "{outdir}/annotation_clones/cluster_clone_counts/{donType}/cluster_clone_counts.ipynb"
    params:
        outdir = lambda wildcards, output: dirname(output.note),
        script = get_lineage_clone_counts_script
    shell: "papermill -p outdir {params.outdir} -p se_cells_meta_f {input.se_meta}  {params.script} {output.note}"


rule plotMarkers:
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
    output:
        note = "{outdir}/annotation_clones/markers/markers.ipynb",
    params:
        rscript = join(ROOT_DIR, "workflow/notebooks/lineage_clones/umap_immune_markers.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
        markers_f = config["markers_f"]
    shell: "papermill -p se_f {input.se_f} -p outdir {params.outdir} -p markers_f {params.markers_f} {params.rscript} {output.note}"

rule overlayClones_umap:
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
    output:
        note = "{outdir}/annotation_clones/umap_clones_overlay/output.ipynb",
    params:
        outdir = lambda wildcards, output: dirname(output.note),
        rscript = join(ROOT_DIR, "workflow/notebooks/lineage_clones/overlayClones_on_umap.ipynb"),
    shell: "papermill -p se_f {input.se_f} -p outdir {params.outdir} {params.rscript} {output.note}"


rule counts_clones:
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
    output:
        multiext("{outdir}/annotation_clones/clone_counts/minCellConds_{min_cells}/",
            "clone_sizes.ipynb", "clone_sizes.csv", "clone_sizes_norm.csv")
    params:
        rscript = join(ROOT_DIR, "workflow/notebooks/lineage_clones/counts_inClones.ipynb"),
        outdir = lambda wildcards, output: dirname(output[0]),
        sample_names = ",".join(config["samples"].index),
    shell:
        "papermill -p se_f {input.se_f} -p outdir {params.outdir} -p sample_names {params.sample_names} -p minCell {wildcards.min_cells} {params.rscript} {output[0]}"


####################################
## Clones-clusters analysis by size.
####################################
rule dominant_clone_clust_in_clone:
    input:
        se_meta = "{outdir}/annotation_clones/se_cells_meta.tsv",
    output:
        note= "{outdir}/annotation_clones/dominant_clone_clust/dominant.ipynb",
        results = report(multiext("{outdir}/annotation_clones/dominant_clone_clust/",
            "dominant_cluster.csv", "cluster_clone_meta.csv", "cluster_clone_counts_normalized.csv",
            "combinedConditions_dominant_cluster.csv", "combinedConditions_cluster_clone_meta.csv", "combinedConditions_cluster_clone_counts_normalized.csv"))
    params:
        outdir = lambda wildcards, output: dirname(output.note),
        script = join(ROOT_DIR, "workflow/notebooks/lineage_clones/clone_cluster_count_embed", "dominant_clust_in_clone.ipynb"),
        samples = ",".join(config["samples"].index),
        cluster_names_f= config.get("cluster_names_f", "None")
    shell:
        "papermill -p outdir {params.outdir} -p se_cells_meta_f {input.se_meta} -p sample_names {params.samples} -p cluster_names_f {params.cluster_names_f:q} {params.script} {output.note}"

rule count_clust_embed:
    input:
        dominant = "{outdir}/annotation_clones/dominant_clone_clust/dominant.ipynb"
    output:
        note = "{outdir}/annotation_clones/clone_clust_embed/tsne_perp{perp}_donperp{donperp}/embed.ipynb",
        sepDonors = "{outdir}/annotation_clones/clone_clust_embed/tsne_perp{perp}_donperp{donperp}/combinedConditions_combinedDonors/umap.pdf",
        combinedDonors = "{outdir}/annotation_clones/clone_clust_embed/tsne_perp{perp}_donperp{donperp}/combinedDonors/umap.pdf"
    #"{outdir}/annotation_clones/clone_clust_umap/"
    params:
        indir = lambda wildcards, input: dirname(input.dominant),
        outdir = lambda wildcards, output: dirname(output.note),
        samples = ",".join(config["samples"].index),
        script = join(ROOT_DIR, "workflow/notebooks/lineage_clones/clone_cluster_count_embed", "run_count_clust_embed.ipynb"),
    shell:  "papermill -p indir {params.indir} -p outdir {params.outdir} -p sample_names {params.samples} -p perplexity {wildcards.perp} {params.script} {output.note}"


rule cluster_clone_hypergeom:
    """Hypergeometric distribution for clones and clusters"""
    input:
        se_meta = "{outdir}/annotation_clones/se_cells_meta.tsv",
    output:
        result="{outdir}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/hypergeom.csv",
        note="{outdir}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/hypergeom.ipynb"
    params:
        rscript = join(ROOT_DIR, "workflow/notebooks/lineage_clones/clones_clusters_hypergeometric.ipynb")
    shell: "papermill -p out_f {output.result} -p se_cells_meta_f {input.se_meta} -p min_clone_size {wildcards.hyperMinCl} -p conds_sep {wildcards.bothConds} -p p_thresh {wildcards.pthresh} {params.rscript} {output.note}"


rule cluster_clone_input_hypergeom:
    """Hypergeometric distribution for clones and clusters"""
    input:
        se_meta = "{outdir}/annotation_clones/se_cells_meta.tsv",
    output:
        result="{outdir}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/input_hypergeom.csv",
        #noIn="{outdir}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/noInput_hypergeom.csv",
        #enrich="{outdir}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/input_hypergeom.csv.png",
        #noIn_enrich="{outdir}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/noInput_hypergeom.csv.png",
        note="{outdir}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/input_hypergeom.ipynb"
    params:
        rscript = join(ROOT_DIR, "workflow/notebooks/lineage_clones/clones_clusters_hypergeometric_sepInput.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
        input_cond = "Input"
    shell: "papermill -p outdir {params.outdir} -p se_cells_meta_f {input.se_meta} -p min_clone_size {wildcards.hyperMinCl} -p conds_sep {wildcards.bothConds} -p p_thresh {wildcards.pthresh} {params.rscript} {output.note}"


def get_hypergeom(wildcards):
    w = wildcards
    hyper = f"{w.outdir}/annotation_clones/hypergeom_clone_clust/mincl.{w.hyperMinCl}_bothConds.{w.bothConds}_p{w.pthresh}/hypergeom.csv"
    if "Input" in config['samples'].index:
        return [hyper, f"{w.outdir}/annotation_clones/hypergeom_clone_clust/mincl.{w.hyperMinCl}_bothConds.{w.bothConds}_p{w.pthresh}/input_hypergeom.csv"]
    return [hyper]


rule run_hypergeom:
    input:
        get_hypergeom
    output:
        "{outdir}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/_complete.txt"
    shell:
        "touch {output}"



rule finalize:
    input:
        dom="{outdir}/annotation_clones/dominant_clone_clust/dominant.ipynb",
        markers="{outdir}/annotation_clones/markers/markers.ipynb",
        cl_sizes=expand("{{outdir}}/annotation_clones/clone_counts/minCellConds_{min_cells}/clone_sizes.ipynb",min_cells=params["min_cells"]),
        cl_lin_sizes = expand("{{outdir}}/annotation_clones/cluster_clone_counts/{donType}/cluster_clone_counts.ipynb",donType=["combinedDonors", "sepDonors"]),
        # hyperg=expand("{{outdir}}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/hypergeom.csv",
        #                hyperMinCl=params["hyperMinCl"], bothConds=params["bothConds"], pthresh=params["p_thresh"]),
        hyperg_f=expand("{{outdir}}/annotation_clones/hypergeom_clone_clust/mincl.{hyperMinCl}_bothConds.{bothConds}_p{pthresh}/_complete.txt",
            hyperMinCl=params["hyperMinCl"], bothConds=params["bothConds"], pthresh=params["p_thresh"]),
        embed = expand("{{outdir}}/annotation_clones/clone_clust_embed/tsne_perp{perp}_donperp{donperp}/embed.ipynb",
                        perp=params_cl_embed["perplexity"], donperp=params_cl_embed["donor_perplexity"]),
        overlay = "{outdir}/annotation_clones/umap_clones_overlay/output.ipynb",
    output:
        "{outdir}/annotation_clones/_nuclear_clones_complete.txt"
    shell: "touch {output}"

