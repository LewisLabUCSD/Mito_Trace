from src.config import ROOT_DIR
import pandas as pd
from os.path import join, dirname
import os
import numpy as np
from icecream import ic


wildcard_constraints:
    variants = "simple|mgatkdonor",
    d = "%d"


samples = config["samples"] #pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)

clones_cfg = config["params"]["clones"]
nclonelist = clones_cfg['vireo']['params']['nclonelist']

########################################################################
## Variant Intermediate step for other workflows.
## Extract pileups for each donor from original filter matrices and run mgatk
########################################################################
# def get_coverage(wildcards):
#     return f"{config['cov_indir']['sample']}/{wildcards.sample}.coverage.txt"



########################################################################
## Vireo A: Directly from multiplex output
########################################################################
rule vireo:
    input:
        "{outdir}/multiplex/multiplex.ipynb" #get_input_multiplex
    output:
        expand("{{outdir}}/clones/variants_simple/vireo/nclones{{nclones}}/donor{d}.labels.png", d=np.arange(config["N_DONORS"])),#, category="lineage"),
        expand("{{outdir}}/clones/variants_simple/vireo/nclones{{nclones}}/donor{d}.variants.labels.png", d=np.arange(config["N_DONORS"])),#, category="lineage"),
        "{outdir}/clones/variants_simple/vireo/nclones{nclones}/cells_meta.tsv"
    params:
        notebook=join("src", "vireo", "2_MT_Lineage_Construct.ipynb"),
        #output_notebook = lambda wildcards, output: join(dirname(output[0]), 'clones.ipynb'),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["N_DONORS"],
        nclones= lambda wildcards: wildcards.nclones #",".join([str(x) for x in nclonelist])
    threads: 8
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p nclones {params.nclones} {params.notebook} {params.OUTDIR}/output.ipynb"


########################################################################
## Vireo B: After multiplexing, separate by donor, grab variants from filters
## that overlap with both conditions, and then run mgatk to call variants again.
##
# Extract pileups for each donor from original filter matrices
########################################################################
rule vireo_mgatkdonor:
    input:
        mtx= "{outdir}/clones/variants_mgatkdonor/vireo/donor{d}/mgatk_donor/cellSNP.tag.AD.mtx",
        cells="{outdir}/clones/variants_mgatkdonor/vireo/donor{d}/mgatk_donor/cells_meta.tsv",
    #output: "clones/variants_mgatkdonor/vireo/clones.ipynb"
    output:
        #"{{outdir}}/clones/variants_mgatkdonor/vireo/nclones{nclones}/clones.ipynb",
        "{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}.labels.png",#, category="lineage", subcategory="donorMGATK"),
        "{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}.variants.labels.png",# category="lineage", subcategory="donorMGATK"),
        "{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}_cells_meta.tsv",
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        #N_DONORS =config["N_DONORS"], #config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "2b_MT_Lineage_Construct_mgatkDonors.ipynb"),
        d = lambda wildcards: wildcards.d,
        #workdir = os.getcwd(),
        nclones= lambda wildcards: wildcards.nclones
    threads: 8
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p donor {params.d} -p nclones {params.nclones} {params.notebook} {params.OUTDIR}/output.ipynb"


rule vireo_mgatkdonor_concat:
    input:
        expand("{{outdir}}/clones/variants_mgatkdonor/vireo/nclones{{nclones}}/donor{d}_cells_meta.tsv",
               d = np.arange(config["N_DONORS"])),
    output:
         "{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/cells_meta.tsv",
    run:
        all = []
        for i in input:
            all.append(pd.read_csv(i, sep='\t'))
        all = pd.concat(all, ignore_index=True).sort_values(["donor", "lineage", "donor_index", "lineage_index"])
        if 'level_0' in all.columns:
                all = all.drop('level_0', axis=1)
        ic('all')
        ic(all.head())
        all.to_csv(output[0], sep='\t', index=False)

# ## Enrichment
# rule vireo_enrichment:
#     version: '1.0' # no +1 norm error
#     input:
#         "{outdir}/clones/variants_{variants}/vireo/nclones{nclones}/cells_meta.tsv"
#     output:
#         report("{outdir}/clones/variants_{variants}/vireo/nclones{nclones}/enrichment/volcano_Fisher_foldNorm.png", category="enrichment"),
#     params:
#         #clones_indir = lambda wildcards, input: dirname(input[0]),#lambda wildcards, input: dirname(input[0]),
#         OUTDIR = lambda wildcards, output: dirname(output[0]),
#         script = join("src", "lineage", "lineage_enrichment.py"),
#         samples=",".join(samples.index)
#     shell: "python {params.script} {input} {params.OUTDIR} {params.samples}"


# Bring enrichment results into same space
rule vireo_process:
    input:
        expand("{{outdir}}/clones/variants_{{variants}}/vireo/nclones{nclones}/enrichment/volcano_Fisher_foldNorm.png",
                    nclones=nclonelist)
    output:
        temp("{outdir}/clones/variants_{variants}/vireo/temp/.tmp")

    shell: "touch {output}"

################################################
## Till here
################################################

# Sort of like an all rule to bring downstream results together.
rule complete_lineage:
    input:
        "{outdir}/clones/variants_{variants}/{method}/temp/.tmp",
    output: "{outdir}/clones/variants_{variants}/{method}/.completed"
    shell: "touch {output}"


