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


rule multiplex:
    input: "{outdir}/cellSNP.tag.AD.mtx"
    output:
        note="{outdir}/multiplex/multiplex.ipynb",
        results = report(multiext("{outdir}/multiplex/", "multiplex_AF_SNPs_all_afFilt.png", "multiplex_clusters_all.labels.png"))
    params:
        N_DONORS=config["N_DONORS"],
        notebook=join("src", "vireo", "1_MT_Donors_multiplex.ipynb" ),
        sample_names= ','.join(samples["sample_name"].values), # make it as a list
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        to_elbo = False
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p sample_names {params.sample_names} -p to_elbo {params.to_elbo} {params.notebook} {output.note}"


rule donors_plotAF:
    input:
        #clone="data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/lineages/lineage.ipynb",
        "{outdir}/multiplex/multiplex.ipynb" #"data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/multiplex.ipynb"
    output:
        report(expand("{{outdir}}/multiplex/dendrograms/figures/donor{d}_dendrogram.png",
               d=np.arange(config["N_DONORS"]))),
        note="{{outdir}}/multiplex/dendrograms/multiplex_dendro.ipynb",
    params:
        #INDIR = lambda wildcards, input: dirname(input[0]),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: Path(abspath(output[0])).parents[1], #go to the 'dendrograms' folder
        out_notebook = lambda wildcards, output: join(Path(abspath(output[0])).parents[1], "af_dendro.ipynb"),
        N_DONORS=config['N_DONORS'], #config["multiplex"]["N_DONORS"],
        sample_names= ",".join(config['samples'].index), # make it as a list
        notebook=join("src", "vireo", "3_MT_Donors_Dendrogram.ipynb"),
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} {params.notebook} {output.note}"


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
