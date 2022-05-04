from os.path import join, dirname
from src.config import ROOT_DIR
from snakemake.utils import min_version
min_version("6.0")
# rule mtclone_:
#     shell:


rule mtclones:
    input:
        don_dir = "{outdir}/multiplex/clones_{variant}/",
    output:
          #note = "{outdir}/multiplex/clones_simpleUnion/mt_clones_thresh/out.ipynb",
          outdir  = "{outdir}/multiplex/clones_{variant}/mt_clones_thresh/",
    params:
        script = join(ROOT_DIR, "workflow/notebooks/mt_as_clones/mt_to_clones_informative.py.ipynb"),
        outdir = lambda wildcards, output: dirname(output),
        mtclone_params = config["mtclone"]["params"],
    log:
        note = "{outdir}/multiplex/clones_simpleUnion/mt_clones_thresh/out.ipynb",
    notebook:
        join(ROOT_DIR, "workflow/notebooks/mt_as_clones/mt_to_clones_informative.ipynb")

    #shell: "papermill -p don_dir {input.don_dir} -p params {params.mtclone_params} {params.script} {params.outdir}"


use_small = False
to_plots = False
rule nuclear_and_mtclone_counts:
    input:
        indir = "{outdir}/multiplex/clones_{variant}/mt_clones_thresh",
         # need to change se_meta to annotation input
         # Get the kparam before
        se_meta = "{outdir}/clones/variants_{variant}/knn/kparam_{kparam}/gff_A2_black/annotation_clones/se_cells_meta_labels.tsv"
    output:
        #note = "{outdir}/multiplex/clones_{variant}/anno_mt_clones_thresh/anno_mt_clones.ipynb",
        note = "{outdir}/multiplex/clones_{variant}/anno_mt_clones_thresh/af.{af}_othaf.{othaf}_cov.{cov}_othcov.{othcov}_ncells.{ncells}_othncells.{othncells}_mean.{mean}/anno_mt_clones.ipynb",
    params:
        script = join(ROOT_DIR, "workflow/notebooks/mt_as_clones/nuclear_and_mt_to_clones_informative.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
        N_DONORS= config["N_DONORS"],
    shell: """ \
           papermill -p indir {input.indir} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} \
           -p af_t {wildcards.af} -p oth_af_t {wildcards.othaf} -p cov_t {wildcards.cov} -p oth_cov_t {wildcards.othcov} \
           -p ncells {wildcards.ncells} -p oth_ncells {wildcards.oth_ncells} -p mean_pos_cov {wildcards.mean_pos_cov} {params.script} {output.note}"""



rule nuclear_mtclones_af:
    input:
        indir = "{outdir}/multiplex/clones_{variant}/",
        se_meta = "{outdir}/clones/variants_{variant}/knn/kparam_{kparam}/gff_A2_black/annotation_clones/se_cells_meta_labels.tsv"
    output:
        note = "{outdir}/clones/variants_{variants}/knn/kparam_{kparam}/gff_A2_black/mt_as_clones/out.ipynb"
    params:
        outdir = lambda wildcards, output: dirname(output.note),
        script = join(ROOT_DIR,"workflow/notebooks/mt_as_clones/mt_nuclear_af.ipynb")
    shell: "papermill -p indir {input.indir} -p se_meta {input.se_meta} -p outdir {params.outdir} {params.script} {output.note}"


rule finalize:
    output:
        "{outdir}/clones/variants_{variants}/mt_as_clones/.complete.txt"
    shell: "touch {output}"