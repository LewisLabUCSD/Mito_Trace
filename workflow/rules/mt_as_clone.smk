from os.path import join, dirname
from src.config import ROOT_DIR
import numpy as np
from snakemake.utils import min_version
min_version("6.0")
# rule mtclone_:
#     shell:

def get_se_meta():
    return join(config["anno_res"], f"gff_{config['gff']}", "mergedSamples","se_cells_meta.tsv")


rule mtclones_informative:
    input:
        #don_dir = "{outdir}/multiplex/clones_{variants}/
        don_dir = expand("{{outdir}}/multiplex/clones_{{variants}}/donor{d}/af.tsv",d=np.arange(config["N_DONORS"]))
    output:
          #note = "{outdir}/multiplex/clones_simpleUnion/mt_clones_thresh/out.ipynb",
          out = "{outdir}/mt_as_clones/variants_{variants}/mt_clones_thresh/.complete",
    params:
        outdir = lambda wildcards, output: dirname(output.out),
        mtclone_params = config["mt_as_clones"]["params"],
    log:
        notebook = "{outdir}/logs/mt_as_clones/variants_{variants}/mt_clones_thresh/out.ipynb",
    notebook:
        join(ROOT_DIR, "workflow/notebooks/mt_as_clones/mt_to_clones_informative.py.ipynb")

    #shell: "papermill -p don_dir {input.don_dir} -p params {params.mtclone_params} {params.script} {params.outdir}"

rule plot_mtclones_informative:
    input:
        #don_dir = "{outdir}/multiplex/clones_{variants}/
        indir = "{outdir}/mt_as_clones/variants_{variants}/mt_clones_thresh/.complete",
    output:
          note = "{outdir}/mt_as_clones/variants_{variants}/mt_clones_thresh/params_plot.ipynb",
          #fig = expand("{{outdir}}/mt_as_clones/variants_{{variants}}/mt_clones_thresh/donor{d}_params_scatter.png",d=np.arange(config["N_DONORS"]))
    params:
        indir = lambda wildcards, input: dirname(input.indir),
        outdir = lambda wildcards, output: dirname(output.note),
        script = join(ROOT_DIR, "workflow/notebooks/mt_as_clones/plot_mt_to_clones_informative.py.ipynb"),
        N_DONORS = config["N_DONORS"]
    log:
        notebook = "{outdir}/logs/mt_as_clones/variants_{variants}/mt_clones_thresh/out.ipynb",
    shell: "papermill -p indir {params.indir} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} {params.script} {output.note}"


use_small = False
to_plots = False
best_p = config["mt_as_clones"]["best_params"]
print("n donors", config["N_DONORS"])
rule nuclear_and_mtclone_counts:
    input:
        indir = "{outdir}/mt_as_clones/variants_{variants}/mt_clones_thresh/.complete", # "{outdir}/multiplex/clones_{variants}/mt_clones_thresh",
         # need to change se_meta to annotation input
         # Get the kparam before
        se_meta = expand("{{outdir}}/anno_multiplex/gff_{gff}/se_cells_meta_labels.tsv", gff=[config["gff"]]) #get_se_meta() #"{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/gff_A2_black/annotation_clones/se_cells_meta_labels.tsv"
    output:
        #note = "{outdir}/multiplex/clones_{variants}/anno_mt_clones_thresh/anno_mt_clones.ipynb",
        note = "{outdir}/mt_as_clones/variants_{variants}/lineage_clusters/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/anno_mt_clones.ipynb",
    params:
        indir  = lambda wildcards, input: dirname(input.indir),
        outdir = lambda wildcards, output: dirname(output.note),
        script = join(ROOT_DIR, "workflow/notebooks/mt_as_clones/nuclear_and_mt_to_clones_informative.ipynb"),
        N_DONORS= config["N_DONORS"],
    #shell: "papermill -p indir {input.indir} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} -p best_p {params.script} {output.note}"
    shell: """ \
           papermill -p indir {params.indir} -p se_meta {input.se_meta} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} \
           -p af_t {wildcards.af} -p oth_af_t {wildcards.othaf} -p cov_t {wildcards.cov} -p oth_cov_t {wildcards.othcov} \
           -p ncells {wildcards.ncells} -p oth_ncells {wildcards.othncells} -p mean_pos_cov {wildcards.mean} {params.script} {output.note}"""


rule best_params_out:
    input:
        indir = "{outdir}/mt_as_clones/variants_{variants}/mt_clones_thresh/.complete", # "{outdir}/multiplex/clones_{variants}/mt_clones_thresh",
    output:
        note = "{outdir}/mt_as_clones/variants_{variants}/bestparams/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/best_params_save.ipynb",
    params:
        indir  = lambda wildcards, input: dirname(input.indir),
        outdir = lambda wildcards, output: dirname(output.note),
        script = join(ROOT_DIR, "workflow/notebooks/mt_as_clones/best_params_mt_to_clones_informative.ipynb"),
        N_DONORS= config["N_DONORS"],
    shell: """ \
           papermill -p indir {params.indir} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} \
           -p af_t {wildcards.af} -p oth_af_t {wildcards.othaf} -p cov_t {wildcards.cov} -p oth_cov_t {wildcards.othcov} \
           -p ncells {wildcards.ncells} -p oth_ncells {wildcards.othncells} -p mean_pos_cov {wildcards.mean} {params.script} {output.note}"""


rule mtclone_bin_dendro:
    input:
        #indir = "{outdir}/mt_as_clones/variants_{variants}/mt_clones_thresh/.complete", # "{outdir}/multiplex/clones_{variants}/mt_clones_thresh",
        indir = "{outdir}/mt_as_clones/variants_{variants}/bestparams/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/best_params_save.ipynb",
    output:
        note = "{outdir}/mt_as_clones/variants_{variants}/dendro/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/dendro_mt_clones.ipynb",
    params:
        indir  = lambda wildcards, input: dirname(input.indir),
        outdir = lambda wildcards, output: dirname(output.note),
        script = join(ROOT_DIR, "workflow/notebooks/mt_as_clones/dendro_mt_to_clones_informative.ipynb"),
        N_DONORS= config["N_DONORS"],
    #shell: "papermill -p indir {input.indir} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} -p best_p {params.script} {output.note}"
    shell: """ papermill -p indir {params.indir} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} {params.script} {output.note}"""

    # shell: """ \
    #        papermill -p indir {params.indir} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} \
    #        -p af_t {wildcards.af} -p oth_af_t {wildcards.othaf} -p cov_t {wildcards.cov} -p oth_cov_t {wildcards.othcov} \
    #        -p ncells {wildcards.ncells} -p oth_ncells {wildcards.othncells} -p mean_pos_cov {wildcards.mean} {params.script} {output.note}"""


rule nuclear_mtclones_af:
    input:
        indir = expand("{{outdir}}/multiplex/clones_{{variants}}/donor{d}/af.tsv",d=np.arange(config["N_DONORS"])), #"{outdir}/multiplex/clones_{variants}/",
        se_meta = expand("{{outdir}}/anno_multiplex/gff_{gff}/se_cells_meta_labels.tsv", gff=config["gff"]) #"{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/gff_A2_black/annotation_clones/se_cells_meta_labels.tsv"
    output:
        note = "{outdir}/mt_as_clones/variants_{variants}/anno_mt_af/out.ipynb"
    params:
        don_dir = lambda wildcards, input: dirname(dirname(input.indir[0])),
        outdir = lambda wildcards, output: dirname(output.note),
        script = join(ROOT_DIR,"workflow/notebooks/mt_as_clones/mt_nuclear_af.ipynb"),
    shell: "papermill -p indir {params.don_dir} -p se_meta {input.se_meta} -p outdir {params.outdir} {params.script} {output.note}"



rule finalize:
    input:
        "{outdir}/mt_as_clones/variants_{variants}/mt_clones_thresh/.complete",
        "{outdir}/mt_as_clones/variants_{variants}/anno_mt_af/out.ipynb",
        "{outdir}/mt_as_clones/variants_{variants}/mt_clones_thresh/params_plot.ipynb",
        expand("{{outdir}}/mt_as_clones/variants_{{variants}}/lineage_clusters/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/anno_mt_clones.ipynb",
                af= best_p["af"], othaf= best_p["oth_af"],
                cov= best_p["cov"], othcov= best_p["oth_cov"],
                ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
                mean= best_p["mean_pos_cov"]),
        expand("{{outdir}}/mt_as_clones/variants_{{variants}}/dendro/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/dendro_mt_clones.ipynb",
                af= best_p["af"],othaf= best_p["oth_af"],
                cov= best_p["cov"], othcov= best_p["oth_cov"],
                ncells= best_p["num_cells"], othncells= best_p["oth_num_cells"],
                mean= best_p["mean_pos_cov"]),
        #"{outdir}/mt_as_clones/variants_{variants}/anno_mt_clones_thresh/manual_params/anno_mt_clones.ipynb",
    output:
        out = "{outdir}/mt_as_clones/variants_{variants}/.complete.txt"
    shell: "touch {output}"