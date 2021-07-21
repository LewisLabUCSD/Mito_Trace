import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from os.path import dirname, join
from src.utils.parse_config import read_config_file
import matplotlib as mpl
mpl.use('Agg')
from sklearn.metrics import silhouette_samples

def get_all_in(wc):
    return

def get_run(wc):
    return allConfig[wc.run]["mtscATAC_OUTDIR"]

allConfig = {}
for c in config["configfiles"]:
    allConfig[c] = read_config_file(config["configfiles"][c])


rule all:
    input:
        expand("{results}/10x_out/10x_qcAll.png",
               results=config["results"])


rule cellr_in_ind:
    #input: get_run
    output:
        report("{results}/10x_out/10x_qcAll.png"),
        # expand("{{results}}/10x_out/10 x_qc_{run}.png",
        #        run=list(config["configfiles"].keys()))
    params:
        outdir= lambda wildcards, output: dirname(output[0])
    #run:


rule cells_mt_mean_and_cv:
    input:
        "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt"
    output:
        report("{results}/countsFilters/cells_dist.png"),
        report("{results}/countsFilters/positions_dist.png")


rule mgatk_cells_mt_mean_and_cv:
    input:
        "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.variant.rds"
    output:
        report("{results}/mgatk_cells/cells_dist.png"),
        report("{results}/mgatk_cells/vars_dist.png")


rule cluster_results:
    input:
        "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/cellSNP.tag.AD.mtx",
        "{results}/data/merged/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/multiplex/multiplex.ipynb"
    output:
        report("{results}/clusters/silhouette.png"),
        report("{results}/clusters/tree_dist.png"),
        report("{results}/clusters/top_transitions.png"),
    params:
        label_col = "donor"
    run:
        s_all = {}
        curr=""
        labels = pd.read_csv(input[0], sep="\t")
        df = pd.read_csv(input[1])
        s = silhouette_samples(df, labels[params.label_col].values, metric='euclidean', sample_size=None, random_state=None)
        s_all[curr] = pd.DataFrame(s,columns=["Score"]).assign(exp=curr)
        sns.violinplot(s)
        #sns.violinplot(pd.concat(s_all), hue="exp")

# rule cellr_input_all:
#     input:
#         get_all_in
#     output:
#         report("{results}/10x_out/10x_qc.png"),
#     params:
#         outdir= lambda wildcards, output: dirname(output[0])
#     run:
#         all_samples = pd.concat((samples, samples_Lar))
#         all_samples["data"] = "Lareau et al 2020"
#         all_samples.loc[samples.index,"data"]="Current"
#
#         all_samples_long = all_samples.reset_index().melt(id_vars=["sample", "data"], value_vars=cols)
#         all_samples_long = all_samples_long.dropna(axis=0)
#         all_samples_long["data"] = all_samples_long["data"].astype(str)
#
#
#         sns.catplot(x="sample", y="value", col="variable", hue="data",
#                    data=all_samples_long, kind='bar', col_wrap=3,
#            sharey=False, palette="colorblind")
#         plt.savefig(join(outdir, "QC_lareau.png"))
