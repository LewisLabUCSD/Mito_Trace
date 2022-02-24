from src.config import ROOT_DIR
import pandas as pd
from os.path import join, dirname
import numpy as np
from icecream import ic


# wildcard_constraints:
#     d = "%d"


# samples = config["samples"] #pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)
#
# clones_cfg = config["params"]["clones"]
# nclonelist = clones_cfg['vireo']['params']['nclonelist']

########################################################################
## Variants C: After multiplexing, separate by donor, keep variants from filters
## that overlap with both conditions, and then use those for clonal calling (may include a lot of variants).
##
# Extract pileups for each donor from original filter matrices
########################################################################
rule knn_simpleUnion:
    input:
        af ="{outdir}/multiplex/clones_simpleUnion/donor{d}/af.tsv"
    output:
        cells = "{outdir}/clones/variants_simpleUnion/knn/kparam_{kparam}/donor{d}/cells_meta.tsv",
        fig = "{outdir}/clones/variants_simpleUnion/knn/kparam_{kparam}/donor{d}/donor{d}.variants.labels.png"
    # fig=report("{outdir}/clones/variants_mgatkdonor/knn/kparam_{kparam}/donor{d}/donor{d}.variants.labels.png",
    #        category="enrichment")
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        name = lambda wildcards: f"donor{wildcards.d}", #"/data2/mito_lineage/data/processed/mttrace/TcellDupi_may17_2021/MTblacklist/pre/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/pre.variant.rds"
        donor = lambda wildcards: wildcards.d,
        outdir = lambda wildcards, output: dirname(output[0]),
        kparam = lambda wildcards: wildcards.kparam,
        note = join(ROOT_DIR, "R_scripts", "knn_clones_init.ipynb"),
    shell:
        "papermill -k ir -p indir {params.indir} -p name {params.name} -p donor {params.donor} -p outdir {params.outdir} -p kparam {params.kparam} {params.note} {params.outdir}/output.ipynb"

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


rule knn_simpleUnion_concat:
    input:
        cells = expand("{{outdir}}/clones/variants_simpleUnion/knn/kparam_{{kparam}}/donor{d}/cells_meta.tsv",
                        d=np.arange(config["N_DONORS"])),
    output:
        cells_all = "{outdir}/clones/variants_simpleUnion/knn/kparam_{kparam}/cells_meta.tsv"
    run: concat_knn(input.cells, output[0])


rule vireo_simpleUnion:
    input:
        af ="{outdir}/multiplex/clones_simpleUnion/donor{d}/af.tsv"
    output:
        #"{{output}}/clones/variants_mgatkdonor/vireo/nclones{nclones}/clones.ipynb",
        "{outdir}/clones/variants_simpleUnion/vireo/nclones{nclones}/donor{d}.labels.png",#, category="lineage", subcategory="donorMGATK"),
        "{outdir}/clones/variants_simpleUnion/vireo/nclones{nclones}/donor{d}.variants.labels.png",# category="lineage", subcategory="donorMGATK"),
        "{outdir}/clones/variants_simpleUnion/vireo/nclones{nclones}/donor{d}_cells_meta.tsv",
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        #N_barcodes =config["N_barcodes"], #config["multiplex"]["N_barcodes"],
        notebook=join("src", "vireo", "2b_MT_Lineage_Construct_mgatkbarcodes.ipynb"),
        d = lambda wildcards: wildcards.d,
        #workdir = os.getcwd(),
        nclones= lambda wildcards: wildcards.nclones
    threads: 8
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p donor {params.d} -p nclones {params.nclones} {params.notebook} {params.OUTDIR}/output.ipynb"


rule vireo_simpleUnion_concat:
    input:
        expand("{{output}}/clones/variants_simpleUnion/vireo/nclones{{nclones}}/donor{d}_cells_meta.tsv",
            d = np.arange(config["N_DONORS"])),
    output:
        "{outdir}/clones/variants_simpleUnion/vireo/nclones{nclones}/cells_meta.tsv",
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