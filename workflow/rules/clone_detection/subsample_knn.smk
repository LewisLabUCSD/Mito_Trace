from src.config import ROOT_DIR
import pandas as pd
from os.path import join, dirname
import os
import numpy as np
from icecream import ic
#np.random.seed(42)
from sklearn.metrics import adjusted_mutual_info_score
from sklearn.metrics import silhouette_score

#samples = config["samples"] #pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)
params = config["clones"]["params"]["subsample"]

# rule all:
#     input:
#         expand("{outdir}/subsamples/pct_{pct}/{i}/donor{d}/knn/kparam_{kparam}/anmi.txt",
#                pct=params["pct"],
#                i=np.arange(params["n_subsamples"]),
#                d=np.arange(config["N_DONORS"])),
#         expand("{outdir}/subsamples/pct_{pct}/{i}/donor{d}/knn/kparam_{kparam}/anmi.txt",
#                pct=params["pct"],
#                i=np.arange(params["n_subsamples"]),
#                d=np.arange(config["N_DONORS"]))

#multiplex/clones_prefilterMerge_impute/donor0/af.tsv
rule subsample_af:
    input:
        af_full = "af.tsv"
    output:
        af_sub = ("{outdir}/subsamples/pct_{pct}/{i}/donor{d}/af.tsv")
    run:
        df = pd.read_csv(input.af_full, sep="\t", index_col=0)
        # keep condittions balanced?
        df.index.name = None
        print(float(wildcards["pct"])/100)
        df.sample(frac=float(wildcards["pct"])/100).to_csv(output.af_sub, sep="\t", index=True)

def get_knn_script(cfg):
    variants_type = cfg.get("variants_type", None)
    if variants_type == "init":
        return join(ROOT_DIR, "R_scripts", "knn_clones_init.ipynb")
    elif variants_type == "mgatkdonor":
        return join(ROOT_DIR, "R_scripts", "knn_clones.ipynb"),
    else:
        return join(ROOT_DIR, "R_scripts", "knn_clones_init.ipynb")

rule knn:
    input:
        af ="{outdir}/subsamples/pct_{pct}/{i}/donor{d}/af.tsv"
    output:
        cells= "{outdir}/subsamples/pct_{pct}/{i}/donor{d}/knn/kparam_{kparam}/cells_meta.tsv",
        fig = "{outdir}/subsamples/pct_{pct}/{i}/donor{d}/knn/kparam_{kparam}/donor{d}.variants.labels.png"
    params:
        indir = lambda wildcards, input: dirname(input.af),
        name = lambda wildcards: f"donor{wildcards.d}",
        outdir = lambda wildcards, output: dirname(output[0]),
        note = join(ROOT_DIR, "R_scripts", "knn_clones_init.ipynb") #get_knn_script(config) #join(ROOT_DIR, "R_scripts", "knn_clones_init.ipynb"),
    priority: 8
    shell:
        "papermill -k ir -p indir {params.indir} -p name {params.name} -p donor {wildcards.d} -p outdir {params.outdir} -p kparam {wildcards.kparam} {params.note} {params.outdir}/output.ipynb"


rule clone_metrics:
    input:
        af ="{outdir}/subsamples/pct_{pct}/{i}/donor{d}/af.tsv",
        cells= "{outdir}/subsamples/pct_{pct}/{i}/donor{d}/knn/kparam_{kparam}/cells_meta.tsv",
    output:
        silhouette = "{outdir}/subsamples/pct_{pct}/{i}/donor{d}/knn/kparam_{kparam}/silhouette.txt",
    run:

        clones = pd.read_csv(input.cells, sep="\t")
        df = pd.read_csv(input.af, sep="\t", index_col=0)

        df = df.loc[clones["ID"]]
        silh = silhouette_score(np.sqrt(df), labels=clones["lineage"])
        with open(output.silhouette, 'w') as f:
            f.write(str(silh))

rule compare_orig:
    input:
        clones_orig = "{outdir}/knn/kparam_{kparam}/cells_meta.tsv",
        clones_sub = "{outdir}/subsamples/pct_{pct}/{i}/donor{d}/knn/kparam_{kparam}/cells_meta.tsv",
    output:
        anmi = "{outdir}/subsamples/pct_{pct}/{i}/donor{d}/knn/kparam_{kparam}/anmi.txt",
    priority: 9
    run:
        clones_orig = pd.read_csv(input.clones_orig, sep="\t")
        clones_orig["donor"] = clones_orig["donor"].astype(str)
        print(clones_orig.head())
        clones_orig = clones_orig[clones_orig["donor"] == str(wildcards.d)]
        clones_sub = pd.read_csv(input.clones_sub, sep="\t")
        inds = list(set(clones_sub["ID"].values).intersection(set(clones_orig["ID"].values)))
        #print('inds', inds[:2])
        clones_orig = clones_orig.set_index("ID").loc[inds,"lineage"].values
        clones_sub  = clones_sub.set_index("ID").loc[inds, "lineage"].values

        anmi = adjusted_mutual_info_score(clones_orig, clones_sub, average_method='arithmetic')
        #print('anmi', anmi)
        with open(output.anmi, 'w') as f:
            f.write(str(anmi))


rule subsamp_concat:
    input:
        silhouette = "{outdir}/subsamples/pct_{pct}/{i}/donor{d}/knn/kparam_{kparam}/silhouette.txt",
        anmi = "{outdir}/subsamples/pct_{pct}/{i}/donor{d}/knn/kparam_{kparam}/anmi.txt",
    output:
        "{outdir}/subsamples/pct_{pct}/{i}/donor{d}/knn/kparam_{kparam}/subsamp.csv",
    priority: 10
    run:
        with open(input.silhouette, 'r') as f:
            silh = f.read()
        with open(input.anmi, 'r') as f:
            anmi = f.read()
        pd.Series([silh, anmi, wildcards.pct], index=["silhouette score", "adjusted nmi", "percent subbsample"]).to_csv(output[0])


rule aggregate:
    input:
        expand("{{outdir}}/subsamples/pct_{pct}/{i}/donor{{d}}/knn/kparam_{{kparam}}/subsamp.csv",
               pct=params["pct"],
               i=np.arange(params["n_subsamples"]),
               )
    output:
        "{outdir}/subsamples/knn/kparam_{kparam}/donor{d}/subsample_metrics.tsv",
    run:
        all = []
        for i in input:
            all.append(pd.read_csv(i, index_col=0).transpose())
        pd.concat(all, axis=0, ignore_index=True).to_csv(output[0], sep="\t")


# rule all_aggregate:
#     input:
#         expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
#                 variants=[x for x in config["params"]["variants"] if x != "simple"], kparam=params_clones["knn"]["params"]["resolution"]),


#rule plots:
    #input:
        #groups = ["subsample_metrics.tsv"]
    #output:
        #out_fig = {outdir}/subsamples/roc_plots.png,
        #out_metrics = {outdir}/subsamples/roc_plots.csv



# def concat_knn(in_files, out_name):
#     all = []
#     cols = ["donor", "lineage"]
#     for i in in_files:
#         #print(i)
#         all.append(pd.read_csv(i,sep='\t'))
#         #print(all[-1].head())
#         if "donor_index" in all[-1].columns.values:
#             cols.append("donor_index")
#         if "lineage_index" in all[-1].columns.values:
#             cols.append("lineage_index")
#     all = pd.concat(all,ignore_index=True).sort_values(cols)
#     if 'level_0' in all.columns:
#         all = all.drop('level_0',axis=1)
#     all.to_csv(out_name,sep='\t',index=False)
#
# rule knn_concat:
#     input:
#         cells= expand("{{outdir}}/knn/kparam_{{kparam}}/donors/donor{d}/cells_meta.tsv",
#                        d=np.arange(config["N_DONORS"])),
#     output:
#         cells_all="{outdir}/knn/kparam_{kparam}/cells_meta.tsv"
#     run: concat_knn(input, output[0])


