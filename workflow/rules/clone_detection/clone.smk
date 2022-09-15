from src.config import ROOT_DIR
import pandas as pd
from os.path import join, dirname
import os
import numpy as np
from icecream import ic


if type(config["dendro_thresh"]) != list:
    dendro_thresh = [config["dendro_thresh"]]
else:
    dendro_thresh = config["dendro_thresh"]


########################################################################
## Variant Intermediate step for other workflows.
## Extract pileups for each donor from original filter matrices and run mgatk
########################################################################
# def get_coverage(wildcards):
#     return f"{config['cov_indir']['sample']}/{wildcards.sample}.coverage.txt"
rule convert_to_af:
    input:
        cells_meta = "{outdir}/cells_meta.tsv",
        counts_in = "af.tsv"
    output:
        note="{outdir}/sc_af/donor{d}/sc_af.ipynb"
        #af="{outdir}/sc_af/donor{d}/af.tsv",
    params:
        note = join(ROOT_DIR, "workflow", "notebooks", "clone_af_dendrograms", "Convert_to_AF.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
        indir = lambda wildcards, input: dirname(input.cells_meta),
        counts_indir = lambda wildcards, input: dirname(input.counts_in),
        var_type = "init"
    shell:
        "papermill -p INDIR {params.indir} -p COUNT_INDIR {params.counts_indir} -p OUTDIR {params.outdir} -p DONOR {wildcards.d}  -p var_type {params.var_type} {params.note} {output.note}"


rule barcodes_btwnClones:
    input:
        cells_meta = "{outdir}/cells_meta.tsv",
        af_note = "{outdir}/sc_af/donor{d}/sc_af.ipynb",
    output:
        note = "{outdir}/barcodes/btwnClones/donor{d}.ipynb",
        # dendros = report("{outdir}/barcodes/btwnClones/donor{d}.na.clust.max2.AF.png",
        #                   category="lineage", subcategory="Clone Barcodes (conditions separate): Donor {d}"),
        # dendrodp  =  report("{outdir}/barcodes/btwnClones/donor{d}.na.clust.max2.DP.png",
        #     category="lineage", subcategory="Clone Barcodes (conditions separate): Donor {d}"),
        # dendrosCond = report("{outdir}/barcodes/btwnClones/donor{d}.NoCondition.na.clust.max2.AF.png",
        #     category="lineage",subcategory="Clone Barcodes: Donor {d}"),
        # dendrodpCond =  report("{outdir}/barcodes/btwnClones/donor{d}.NoCondition.na.clust.max2.DP.png",
        #     category="lineage", subcategory="Clone Barcodes: Donor {d}"),
    params:
        note = join(ROOT_DIR, "workflow", "notebooks", "clone_af_dendrograms", "MT_btwnClones_Barcode.ipynb"),
        indir = lambda wildcards, input: dirname(input.cells_meta),
        outdir = lambda wildcards, output: dirname(output.note),
    shell: "papermill -p INDIR {params.indir} -p OUTDIR {params.outdir} -p DONOR {wildcards.d}  {params.note} {output.note}"



rule barcodes_btwnClones_dendro:
    input:
        cells_meta = "{outdir}/cells_meta.tsv",
        af_note = "{outdir}/sc_af/donor{d}/sc_af.ipynb",
    output:
        note = "{outdir}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.ipynb",
        mean = "{outdir}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.mean.csv",
        res = report(multiext("{outdir}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.",
                                  "clones_dendro.csv", "dendrogram_pvals.txt",
                                  "dendro.NoCondition.max2.AF.png"),
                         category="lineage")

    params:
        note = join(ROOT_DIR, "workflow", "notebooks", "clone_af_dendrograms", "MT_btwnClones_Barcode_dendro.ipynb"),
        indir = lambda wildcards, input: dirname(input.cells_meta),
        outdir = lambda wildcards, output: dirname(output.note),
        #dendro_thresh = dendro_thresh #0.6
    shell: "papermill -p INDIR {params.indir} -p OUTDIR {params.outdir} -p DONOR {wildcards.d} -p dendroThresh {wildcards.dendro_thresh} {params.note} {output.note}"


rule barcodes_inClones:
    input:
        cells_meta = "{outdir}/cells_meta.tsv",
        af_note = "{outdir}/sc_af/donor{d}/sc_af.ipynb",
    output:
        note = "{outdir}/barcodes/inClones/donor{d}.ipynb",
    params:
        note = join(ROOT_DIR, "workflow/notebooks/clone_af_dendrograms", "MT_inClones_Barcode.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
        indir = lambda wildcards, input: dirname(input.cells_meta),
    shell: "papermill -p INDIR {params.indir} -p OUTDIR {params.outdir} -p DONOR {wildcards.d} {params.note} {output.note}"


# rule distinguishing_vars:
#     """ Code tto find enriched barcodes in clones.
#         Replaced by moduel enriched_clone_variants.smk
#     """
#     input:
#         cells_meta = "{outdir}/cells_meta.tsv",
#         af_note = "{outdir}/sc_af/donor{d}/sc_af.ipynb",
#     output:
#         note = "{outdir}/distinct_variants/donor{d}/output.ipynb",
#     params:
#         note = join(ROOT_DIR, "workflow/notebooks/clone_vars", "distinguishing_vars.ipynb"),
#         outdir = lambda wildcards, output: dirname(output.note),
#         indir = lambda wildcards, input: dirname(input.cells_meta),
#     shell: "papermill -p INDIR {params.indir} -p OUTDIR {params.outdir} -p DONOR {wildcards.d} {params.note} {output.note}"
#

rule finalize:
    input:
        #inClones =  expand("{{outdir}}/barcodes/btwnClones/donor{d}.ipynb", d=np.arange(config["N_DONORS"])),
        btwnClones = expand("{{outdir}}/barcodes/inClones/donor{d}.ipynb", d=np.arange(config["N_DONORS"])),
        dendroClones = expand("{{outdir}}/barcodes/btwnClones_dendro_dt_{dendro_thresh}/donor{d}.mean.csv",
                              d=np.arange(config["N_DONORS"]), dendro_thresh=dendro_thresh),
    output:
        out = "{outdir}/barcodes/_clone_complete.txt",
    shell: "touch {output.out}"
