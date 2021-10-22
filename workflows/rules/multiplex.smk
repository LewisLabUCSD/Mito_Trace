from os.path import join, dirname, abspath
import numpy as np
samples = config["samples"]
from snakemake.utils import min_version
min_version("6.0")
from pathlib import Path

# rule all:
#     input:
#         "multiplex/multiplex.ipynb",
#         expand("dendrograms/figures/donor{d}_dendrogram.png", d=np.arange(config["N_DONORS"])),
#         "multiplex/variants/variants.ipynb",
#

def get_mult_input(wc):
    if config["N_DONORS"] == 1:
        return [f"{wc.outdir}/multiplex/single_multiplex.ipynb",
                  f"{wc.outdir}/multiplex/single_multiplex_AF_SNPs_all_afFilt.png",
                  f"{wc.outdir}/multiplex/single_multiplex_clusters_all.labels.png"]

    else: return [f"{wc.outdir}/multiplex/out_multiplex.ipynb",
                  f"{wc.outdir}/multiplex/out_multiplex_AF_SNPs_all_afFilt.png",
                  f"{wc.outdir}/multiplex/out_multiplex_clusters_all.labels.png"]


rule multiplex:
    input: "{outdir}/cellSNP.tag.AD.mtx"
    output:
        note="{outdir}/multiplex/out_multiplex.ipynb",
        results = multiext("{outdir}/multiplex/", "out_multiplex_AF_SNPs_all_afFilt.png",
                           "out_multiplex_clusters_all.labels.png")#, category="Multiplex")
    params:
        N_DONORS=config["N_DONORS"],
        notebook=join("src", "vireo", "1_MT_Donors_multiplex.ipynb" ),
        sample_names= ','.join(samples["sample_name"].values), # make it as a list
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        to_elbo = False
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p sample_names {params.sample_names} -p to_elbo {params.to_elbo} --log-output {params.notebook} {output.note}"


rule multiplex_single:
    """ Copy over the cells meta and make fake figures
    """
    input: "{outdir}/cellSNP.tag.AD.mtx"
    output:
        "{outdir}/multiplex/single_multiplex.ipynb",
        results = multiext("{outdir}/multiplex/", "single_multiplex_AF_SNPs_all_afFilt.png",
                           "single_multiplex_clusters_all.labels.png")#, category="Multiplex")
    params:
        sample_names= ','.join(samples["sample_name"].values), # make it as a list
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        cells_meta = pd.read_csv(join(params.INDIR, "cells_meta.tsv"),
                                 sep='\t')
        new_cells_meta = cells_meta.copy()
        new_cells_meta["donor"] = 0
        new_cells_meta["donor_index"] = np.arange(1, new_cells_meta.shape[0]+1)
        new_cells_meta.to_csv(join(params.OUTDIR, "cells_meta.tsv"), sep='\t')
        new_cells_meta.to_csv(join(params.OUTDIR, "donor0.labels.txt"))
        f = plt.figure() # save blank figure
        print(output.results[0])
        plt.savefig(output.results[0])
        plt.savefig(output.results[1])
        cmd = f"touch {output[0]}"
        os.system(cmd)
        cmd = f"cp {params.INDIR}/cellSNP.base.vcf {params.OUTDIR}/donor0.vcf "
        os.system(cmd)
        cmd = f"cp {params.INDIR}/cellSNP.tag.AD.mtx {params.OUTDIR}/donor0.AD.mtx "
        os.system(cmd)
        cmd = f"cp {params.INDIR}/cellSNP.tag.DP.mtx {params.OUTDIR}/donor0.DP.mtx"
        os.system(cmd)

rule agg_multiplex:
    input: get_mult_input
    output:
        "{outdir}/multiplex/multiplex.ipynb",
        "{outdir}/multiplex/out_multiplex_AF_SNPs_all_afFilt.png",
        "{outdir}/multiplex/out_multiplex_clusters_all.labels.png"
    shell: "cp {input[0]} {output[0]} & cp {input[1]} {output[1]} & cp {input[2]} {output[2]}"


rule donors_plotAF:
    input:
        #clone="data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/lineages/lineage.ipynb",
        "{outdir}/multiplex/multiplex.ipynb"
    output:
        report(expand("{{outdir}}/multiplex/dendrograms/figures/donor{d}_dendrogram.{suf}",
               d=np.arange(config["N_DONORS"]), suf=["png", "depth.png", "withHigh.png"]), category="multiplex"),
    params:
        note = lambda wildcards, output: join(dirname(dirname(output[0])), "multiplex_dendro.ipynb"), # "{{outdir}}/multiplex/dendrograms/multiplex_dendro.ipynb",
        #INDIR = lambda wildcards, input: dirname(input[0]),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: Path(abspath(output[0])).parents[1], #go to the 'dendrograms' folder
        out_notebook = lambda wildcards, output: join(Path(abspath(output[0])).parents[1], "af_dendro.ipynb"),
        N_DONORS=config['N_DONORS'], #config["multiplex"]["N_DONORS"],
        sample_names= ",".join(config['samples'].index), # make it as a list
        notebook=join("src", "vireo", "3_MT_Donors_Dendrogram.ipynb"),
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} {params.notebook} {params.note}"



rule donors_type_variants:
    input:
        "{outdir}/multiplex/multiplex.ipynb"
    output: "{outdir}/multiplex/variants/variants.ipynb",
    params:
        notebook=join("src", "vireo", join("5_MT_Donors_variantTypes.ipynb")),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["N_DONORS"],
        sample_names= ",".join(config['samples'].index)
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} {params.notebook} {output} && jupyter nbconvert --to pdf {output}"
