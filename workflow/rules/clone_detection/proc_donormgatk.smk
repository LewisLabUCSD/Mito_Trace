from src.config import ROOT_DIR
import pandas as pd
from os.path import join, dirname
import numpy as np
from icecream import ic


# wildcard_constraints:
#     variants = "simple|mgatkdonor|init",
#     d = "%d"


# samples = config["samples"] #pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)
#
# clones_cfg = config["params"]["clones"]
# nclonelist = clones_cfg['vireo']['params']['nclonelist']

########################################################################
## Variant Intermediate step for mgatkdonor
## Extract pileups for each donor from original filter matrices and run mgatk
########################################################################
rule donor_mgatk_to_vireoIn:
    input:
        "{outdir}/clones/variants_mgatkdonor/donor{d}/mgatk/d{d}.variant.rds"
    output:
        "{outdir}/clones/variants_mgatkdonor/vireo/donor{d}/mgatk/cellSNP.tag.AD.mtx",
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        donor = lambda wildcards: f"d{wildcards.d}"
    shell:
        "python src/mgatk_to_vireo.py {params.indir} {params.outdir} {params.donor}"

rule donor_copy_cells_meta:
    input: "{outdir}/multiplex/multiplex.ipynb"
    output:
        "{outdir}/clones/variants_mgatkdonor/vireo/donor{d}/mgatk/cells_meta.tsv"
    params:
        cells_meta = lambda wildcards, input: join(dirname(input[0]), "cells_meta.tsv")
    shell: "cp {params.cells_meta} {output}"


########################################################################
## Vireo B
########################################################################
rule vireo_mgatkdonor:
    input:
        mtx= "{outdir}/clones/variants_mgatkdonor/vireo/donor{d}/mgatk/cellSNP.tag.AD.mtx",
        cells="{outdir}/clones/variants_mgatkdonor/vireo/donor{d}/mgatk/cells_meta.tsv",
    #output: "clones/variants_mgatkdonor/vireo/clones.ipynb"
    output:
        #"{{output}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_mgatkdonor/vireo/nclones{nclones}/clones.ipynb",
        "{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}.labels.png",#, category="lineage", subcategory="donorMGATK"),
        "{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}.variants.labels.png",# category="lineage", subcategory="donorMGATK"),
        "{outdir}/clones/variants_mgatkdonor/vireo/nclones{nclones}/donor{d}_cells_meta.tsv",
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        #N_DONORS =config["N_DONORS"], #config["multiplex"]["N_DONORS"],
        notebook=join("src", "vireo", "2b_MT_Lineage_Construct_mgatkbarcodes.ipynb"),
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

########################################################################
## KNN 
########################################################################
rule knn_mgatkdonor:
    input:
        in_vars="{outdir}/clones/variants_mgatkdonor/donor{d}/mgatk/d{d}.variant.rds"
    output:
        cells= "{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
        fig = "{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/donor{d}.variants.labels.png"
    # fig=report("{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/donor{d}.variants.labels.png",
    #        category="enrichment")
    params:
        mgatk_in = lambda wildcards, input: input[0].replace(".variant.rds", ".af.tsv"),
        name = lambda wildcards: f"donor{wildcards.d}", #"/data2/mito_lineage/data/processed/mttrace/TcellDupi_may17_2021/MTblacklist/pre/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/pre.variant.rds"
        donor = lambda wildcards: wildcards.d,
        outdir = lambda wildcards, output: dirname(output[0]),
        kparam = lambda wildcards: wildcards.kparam,
        note = join(ROOT_DIR, "R_scripts", "knn_clones.ipynb"),
    shell:
        "papermill -k ir -p mgatk_in {params.mgatk_in} -p name {params.name} -p donor {params.donor} -p outdir {params.outdir} -p kparam {params.kparam} {params.note} {params.outdir}/output.ipynb"


def concat_knn(in_files, out_name):
    all = []
    cols = ["donor", "lineage"]
    for i in in_files:
        #print(i)
        all.append(pd.read_csv(i,sep='\t'))
        #print(all[-1].head())
        if "donor_index" in all[-1].columns.values:
            cols.append("donor_index")
        if "lineage_index" in all[-1].columns.values:
            cols.append("lineage_index")
    all = pd.concat(all,ignore_index=True).sort_values(cols)
    if 'level_0' in all.columns:
        all = all.drop('level_0',axis=1)
    all.to_csv(out_name,sep='\t',index=False)


rule knn_mgatkdonor_concat:
    input:
        donors = expand("{{outdir}}/clones/variants_mgatkdonor/knn/kparam_{{kparam}}/donor{d}/cells_meta.tsv",
                       d = np.arange(config["N_DONORS"])),
    output:
        cells_all = "{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/cells_meta.tsv"
    run:
        all = []
        cols = ["donor", "lineage"]
        for i in input:
            #print(i)
            all.append(pd.read_csv(i, sep='\t'))
            #print(all[-1].head())
            if "donor_index" in all[-1].columns.values:
                cols.append("donor_index")
            if "lineage_index" in all[-1].columns.values:
                cols.append("lineage_index")
        all = pd.concat(all, ignore_index=True).sort_values(cols)
        if 'level_0' in all.columns:
            all = all.drop('level_0', axis=1)
        all.to_csv(output[0], sep='\t', index=False)
