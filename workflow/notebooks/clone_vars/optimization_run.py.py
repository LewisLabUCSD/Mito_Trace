#!/usr/bin/env python
# coding: utf-8

import matplotlib
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
#matplotlib.rcParams['text.usetex'] = True
from icecream import ic

ic.disable()


#indir = snakemake.input.indir
indir = snakemake.params.indir
outdir = snakemake.params.outdir
donor =  int(snakemake.params.donor)
anno_cells_meta_f = snakemake.input.anno_cells_meta_f  #"/data/Mito_Trace/output/pipeline/v02/CHIP_b1/MTBlacklist_A2/data/merged/MT/cellr_True/numread_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/mgatk/vireoIn/clones/variants_init/knn/kparam_30/gff_A2_black/annotation_clones/se_cells_meta_labels.tsv"
# Objective weights. order of the columns
weights =  snakemake.params.weights
objectives_l = snakemake.params.get("objectives_l", 
                                    ["variants_with_clone_norm_by_1_over_nclones_with_variant", 
                                     "max_clone_ncells_over_nclones", "max_clone_ncells_over_ncells", 
                                     "pct_thresh","other_pct_thresh", 
                                     "n_vars", "obj_nclones_more_than_one_unique"])
ncpus = snakemake.params.get('ncpus', 8)
topn = snakemake.params.get("topn", 16)

to_test = snakemake.params.get("to_test", False)


# try:
#     indir = snakemake.input.indir
#     outdir = snakemake.params.outdir
#     donor =  snakemake.params.donor
#     anno_cells_meta_f = snakemake.input.anno_cells_meta_f  #"/data/Mito_Trace/output/pipeline/v02/CHIP_b1/MTBlacklist_A2/data/merged/MT/cellr_True/numread_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/mgatk/vireoIn/clones/variants_init/knn/kparam_30/gff_A2_black/annotation_clones/se_cells_meta_labels.tsv"
#     # Objective weights. order of the columns
#     weights =  snakemake.params.weights
#     objectives_l = snakemake.params.get("objectives_l", 
#                                         ["variants_with_clone_norm_by_1_over_nclones_with_variant", 
#                                          "max_clone_ncells_over_nclones", "max_clone_ncells_over_ncells", 
#                                          "pct_thresh","other_pct_thresh", 
#                                          "n_vars", "obj_nclones_more_than_one_unique"])
#     ncpus = snakemake.params.get(ncpus, 8)
#     topn = snakemake.params.get(topn, 16)

# except:                                       
#     indir = "/data/Mito_Trace/output/pipeline/v02/CHIP_b1/MTBlacklist_A2/data/merged/MT/cellr_True/numread_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/mgatk/vireoIn/clones/variants_init/knn/kparam_30/"
#     outdir = "/data/Mito_Trace/output/pipeline/v02/CHIP_b1/MTBlacklist_A2/data/merged/MT/cellr_True/numread_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/mgatk/vireoIn/clones/variants_init/knn/kparam_30/distinct_variants/donor0/scrap/"
#     donor = 0
#     anno_cells_meta_f = "/data/Mito_Trace/output/pipeline/v02/CHIP_b1/MTBlacklist_A2/data/merged/MT/cellr_True/numread_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/mgatk/vireoIn/clones/variants_init/knn/kparam_30/gff_A2_black/annotation_clones/se_cells_meta_labels.tsv"
#     # pct_thresh = [0.01, 0.1, 0.25, 0.4, 0.5, 0.75, 0.95]
#     # other_pct_thresh = [0.01, 0.1, 0.25, 0.5]
#     # af_thresh = [0, 0.01, 0.1, 0.25, 0.4]

#     # Objective weights. order of the columns
#     weights = [1,0,0,1,-1, 1, 1] #[1,1,1,1,1] #np.ones([len(objectives),])
#     objectives_l = ["variants_with_clone_norm_by_1_over_nclones_with_variant", 
#                     "max_clone_ncells_over_nclones", "max_clone_ncells_over_ncells", 
#                     "pct_thresh","other_pct_thresh", "n_vars", "obj_nclones_more_than_one_unique"] #"nvars"
#     ncpus=8
#     topn=16


import multiprocessing
import pandas as pd
import numpy as np
import random
import src.clone_variants_optim  as optim

from os.path import join, exists, dirname
from glob import glob
import pickle
import mplh.cluster_help as ch
import mplh.fig_utils as fu

import os
import seaborn as sns
from scipy import sparse
from scipy.io import mmread
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
np.set_printoptions(formatter={'float': lambda x: format(x, '.3f')})


objectives = {ind:x for ind,x in enumerate(objectives_l)}
weights = np.array(weights)
param_names = ["pct_thresh","af_thresh", "other_pct_thresh"]
n_params = 3


pct_thresh = np.arange(0.05, 1, 0.05)
other_pct_thresh = np.arange(0.005, 1, 0.05)
af_thresh = np.arange(0.005, 1, 0.05)

params = {"pct_thresh": pct_thresh,
          "other_pct_thresh": other_pct_thresh,
          "af_thresh": af_thresh,}
assert(len(weights)==len(objectives))




if not exists(outdir):
    os.mkdir(outdir)

# ## Load & preprocess:
# - AF df
# - DP df
# - cells_meta with clone labels. need to create name as donor_lineage
# 
# Remove donor variants (>0.9 in 90% of pop)

af_indir = join(indir, "sc_af", f"donor{donor}")

AF_df = pd.read_csv(join(af_indir, "af.tsv"), index_col=0, sep="\t")
DP_df = pd.read_csv(join(af_indir, "dp.tsv"), index_col=0, sep="\t")

print(AF_df.shape)
print(DP_df.shape)
print("Depth")
print(DP_df.head())
AF_df.head()

cells_meta = pd.read_csv(join(indir, "cells_meta.tsv"), sep='\t', index_col="ID")#.sort_values(["donor", "lineage"])
cells_meta["name"] = cells_meta["donor"].astype(str)+"_"+cells_meta["lineage"].astype(str)
curr_labels = cells_meta[cells_meta["donor"]==donor]


conditions = curr_labels["condition"].unique()
anno_cells = pd.read_csv(anno_cells_meta_f, sep="\t", index_col=0)


## Get donor inds
donor_inds = AF_df.index[((AF_df>0.9).sum(axis=1)>(0.9*AF_df.shape[1]))]


def rm_high(df, thresh, pct_thresh):
    return df.loc[~(((df>thresh).sum(axis=1)>pct_thresh*df.shape[0]))]

def rm_low(df, thresh, pct_thresh):
    return df.loc[~((df<thresh).sum(axis=1)>(pct_thresh*df.shape[1]))]



# parallel setup
from pandarallel import pandarallel
pandarallel.initialize(nb_workers=ncpus) #, progress_bar=True)
from itertools import product



# There are 3 params used for calling the clone
full_params = list(product(*list(params.values())))
full_params = pd.DataFrame(full_params, columns=params.keys())

print(full_params.shape)
full_params.head()

################################################
# Run pipeline and get objective scores
################################################
if to_test:
    results_df = full_params.sample(32).parallel_apply(
        optim.evaluate_series, args=(AF_df, DP_df, curr_labels), axis=1)
else:
    if full_params.shape[0]>=10000:
        results_df = full_params.sample(10000).parallel_apply(optim.evaluate_series, args=(AF_df, DP_df, curr_labels), axis=1)
    else:
        results_df = full_params.parallel_apply(optim.evaluate_series, args=(AF_df, DP_df, curr_labels), axis=1)
    #results_df = full_params.sample(100).parallel_apply(evaluate_series, args=(AF_df, DP_df, curr_labels), axis=1)

print((results_df>0).any())
drop_inds = results_df.loc[(results_df==0).all(axis=1)].index
results_norm = optim.set_multi(results_df, weights)
rank_df = optim.set_multi_rank(results_norm, weights)

drop_results = results_norm.loc[results_norm["multi"].isnull()]
results_norm = results_norm.loc[~(results_norm["multi"].isnull())]

########################
## Save the objective results
########################
results_norm.to_csv(join(outdir, "objectives_norm.csv"))
results_df.to_csv(join(outdir, "objectives.csv"))
rank_df.to_csv(join(outdir, "objectives_rank.csv"))
#results_norm.loc[results_norm["multi"] == np.nan]

########################
# Plot distribution results
########################
sns.displot(results_norm["multi"])
plt.title("multiobjective function (want to maximize)")
plt.savefig(join(outdir, "loss_multi.pdf"))

sns.displot(rank_df["multi"])

sns.displot(rank_df["variants_with_clone_norm_by_1_over_nclones_with_variant"])
plt.title("objective: variants_with_clone_norm_by_1_over_nclones_with_variant (want to maximize)")
plt.savefig(join(outdir, "loss_variants_with_clone_norm_by_1_over_nclones_with_variant.pdf"))

sns.displot(rank_df["obj_nclones_more_than_one_unique"])
plt.title("objective: Number of clones with at least 2 unique variants (want to maximize)")
plt.savefig(join(outdir, "loss_variants_with_clone_norm_by_1_over_nclones_with_variant.pdf"))


########################
# Get the top n results
########################
def get_top_n_results(results_df, rank_df, n=12):
    filt_rank = rank_df.sort_values(by=["multi"])[::-1].iloc[:n]
    filt_results = results_df.loc[filt_rank.index]
    return filt_rank, filt_results


full_params.to_csv(join(outdir, "params.csv"))

filt_rank, filt_results = get_top_n_results(results_norm, rank_df, n=topn)
filt_results.columns = [f"{x}_obj" for x in filt_results.columns]
filt_results = pd.merge(filt_results, full_params, left_index=True, right_index=True, how="left")
filt_rank = filt_rank.loc[filt_results.index]

all_df = []
all_objs = {}
for ind, val in filt_results.iterrows():
    print(ind)
    obj_out, data = optim.evaluate_series(val, AF_df, DP_df, curr_labels, return_data=True)
    all_df.append(data["all_unique_df"])
    all_objs[ind] = obj_out 
all_df = pd.concat(all_df)

heatmap_input = all_df[["n_cells", "variant"]].reset_index().pivot(index="id", columns="variant", values="n_cells").fillna(0).astype(int)
meta_df = all_df[["af_thresh", "other_pct_thresh", "pct_thresh", "clone"]]
meta_df = meta_df.loc[~(meta_df.index.duplicated())]
meta_df = meta_df.sort_values(["af_thresh","pct_thresh", "other_pct_thresh", "clone"])
heatmap_input = heatmap_input.loc[meta_df.index]

# Get the variants based on total number of cells across parameters
heatmap_input = heatmap_input.loc[:,heatmap_input.sum().sort_values()[::-1].index]
variants_order = heatmap_input.columns


clone_sums = meta_df.groupby("clone").apply(optim.clone_sum, heatmap_input)
clone_sums = clone_sums.loc[:, clone_sums.sum().sort_values()[::-1].index]
clones_order = clone_sums.index

all_df["params"] = all_df.apply(optim.params_to_str, axis=1, args=(param_names,))
filt_results["params"] = filt_results.apply(optim.params_to_str, axis=1, args=(param_names,))
filt_results["params_multi"] = filt_results.apply(optim.params_and_multi_str, axis=1)
tmp = filt_results.set_index("params")
all_df["multi_obj"] = all_df.apply(lambda x: tmp.loc[x["params"], "multi_obj"], axis=1)
del tmp                               

all_df["params_multi"] = all_df.apply(optim.params_and_multi_str, axis=1)


fg = sns.FacetGrid(data=all_df.reset_index(), height=4, sharey=False, sharex=False,
                   col="params", col_wrap=4, col_order=filt_results["params"].values, margin_titles=True)

fg.map_dataframe(optim.draw_heatmap, 'variant','clone', 'log2_n_cells',
                 clones_order=clones_order, variants_order=variants_order)#, cbar=False)

#fg.set_titles(row_template = 'other_pct_thresh: {row_name}', col_template = 'pct_thresh: {col_name}')
fg.fig.suptitle(f"Best parameter combinations shown in order")
fg.fig.subplots_adjust(top=0.9, hspace = 0.8)

plt.title("multiobjective function (want to maximize)")
#plt.savefig(join(outdir, "top_param_results.pdf"))
plt.savefig(join(outdir, "top_param_results.pdf"), dpi=300)



# best_params = (filt_results.sort_values("multi_obj", ascending=False)).iloc[0]
# best_params = pd.DataFrame(best_params).transpose()
# best_params.index = ["objective_scores"]
# best_params.loc["weight"] = None
# for obj, w in zip(objectives_l, weights):
#     best_params.loc["weight", f"{obj}_obj"] = w
#
# out_df = all_df[all_df['params'] == best_params.loc["objective_scores", "params"]]
#
# clone_var_table = (out_df.pivot(index= 'variant',columns='clone', values='log2_n_cells').fillna(0))
#
#
# clones_keep = clone_var_table.loc[:, ~((clone_var_table==0).all(axis=0))].columns
# vars_keep = clone_var_table.loc[~((clone_var_table==0).all(axis=1))].index
#
#
# sns.clustermap(clone_var_table)
# plt.title(best_params.loc["objective_scores", "params_multi"])
# plt.savefig(join(outdir, "best_params.pdf"))
#
# sns.clustermap(clone_var_table.loc[vars_keep,clones_keep])
# plt.title(best_params.loc["objective_scores", "params_multi"])
# plt.savefig(join(outdir, "best_params_filt.pdf"))
#
#
# # ## Save clone-variant table and the parameters
# clone_var_table.to_csv(join(outdir, "best_params_clone_vars.csv"))
# clone_var_table.loc[vars_keep,clones_keep].to_csv(join(outdir, "best_params_filt_clone_vars.csv"))
# best_params.to_csv(join(outdir, "best_params.csv"))
#
#
# filt_curr_labels = curr_labels[curr_labels["name"].isin(clones_keep)]
# #anno_cells = anno_cells.loc[anno_cells["ID"].isin(filt_curr_labels.index)]
# #out_cells_meta = anno_cells.loc[anno_cells["ID"].isin(filt_curr_labels.index)]
#
# # overlap cells of anno and curr labels
# cells_to_keep = set(anno_cells["ID"].values).intersection(set(filt_curr_labels.index))
# out_cells_meta = anno_cells.loc[anno_cells["ID"].isin(cells_to_keep)]
# #out_cells_meta =  out_cells_meta.reset_index().set_index("ID")
#
# # out_AF_df = AF_df.loc[vars_keep, out_cells_meta.index]
# # out_DP_df = DP_df.loc[vars_keep, out_cells_meta.index]
# print(out_cells_meta.head())
# out_AF_df = AF_df.loc[vars_keep, out_cells_meta["ID"]].transpose()
# out_DP_df = DP_df.loc[vars_keep, out_cells_meta["ID"]].transpose()
#
#
# print(out_cells_meta.shape)
# print(out_AF_df.shape)
# print(out_DP_df.shape)
#
# assert((out_AF_df.index==out_DP_df.index).all())
# assert((out_AF_df.columns==out_DP_df.columns).all())
# assert((out_AF_df.index==out_cells_meta["ID"]).all())
#
#
# # ## save cells-meta, af and dp
# # In[ ]:
# out_cells_meta["ID"]
#
# out_cells_meta.to_csv(join(outdir, "cells_meta.tsv"),sep="\t")
# out_AF_df.to_csv(join(outdir, "af.tsv"), sep="\t")
# out_DP_df.to_csv(join(outdir, "dp.tsv"), sep="\t")



# # Get the clones based on total number of cells across parameters
# def clone_sum(val):
#     #print(val)
#     return(heatmap_input.loc[val.index].sum())
# def params_to_str(ser, param_names):
#     name = ""
#     for p in param_names:
#         name = name + f"{p}={ser[p]:.3f}\n"
#     return name
#
#
# def params_and_multi_str(ser):
#     param_str = ser["params"]
#     name = f"params:\n{param_str.strip()}\nObjective score={ser['multi_obj']} (want to maximize)"
#     return name

#
# def get_unique_variants(cln_af, other_af, pct_thresh, af_thresh, other_pct_thresh):
#     """ gets the distinct variants in a clone.
#     """
#     n_thresh = pct_thresh*cln_af.shape[1]
#     n_oth_thresh = other_pct_thresh*other_af.shape[1]
#     bin_cln = cln_af>af_thresh
#     bin_other = other_af>af_thresh
#     cells_above = bin_cln.sum(axis=1)
#     pct_cells_above = cells_above/bin_cln.shape[1]
#     up_vars = bin_cln.loc[cells_above > n_thresh].index
#     cells_other_above = bin_other.sum(axis=1)
#     pct_cells_other_above = cells_other_above/bin_other.shape[1]
#     up_oth_vars = bin_other.loc[cells_other_above > n_oth_thresh].index
#     uniq_vars = list(set(up_vars) - set(up_oth_vars))
#     out = pd.DataFrame(index=uniq_vars, data={"n_cells":cells_above.loc[uniq_vars].values,
#                                               "n_other_cells": cells_other_above.loc[uniq_vars].values,
#                                               "pct_above": pct_cells_above,
#                                               "pct_other_above": pct_cells_other_above})
#     out["pct_thresh"] = pct_thresh
#     out["af_thresh"] = af_thresh
#     out["other_pct_thresh"] = other_pct_thresh
#     return out
#
#
# def get_clones_unique_variants(solution, data):
#     all_unique_df = []
#     pct_thresh, af_thresh, other_pct_thresh = solution["pct_thresh"], solution["af_thresh"], solution["other_pct_thresh"] #solution[0], solution[1], solution[2]
#     curr_labels = data["curr_labels"]
#     AF_df = data["AF_df"]
#     DP_df = data["DP_df"]
#     for cln, val in curr_labels.groupby("name"):
#         ic(cln)
#         cln_af = AF_df.loc[:, val.index]
#         other_af = AF_df.loc[:, curr_labels.drop(val.index).index]
#         curr_dp = DP_df.loc[:, val.index]
#         curr_labs = curr_labels[curr_labels.index.isin(cln_af.columns)]
#         ic(cln_af.shape)
#         unique_df = get_unique_variants(cln_af, other_af, pct_thresh, af_thresh, other_pct_thresh)
#         unique_df["clone"] = cln
#         unique_df["id"] = unique_df["clone"] + "_" + unique_df["pct_thresh"].astype(str)+ "_" + unique_df["af_thresh"].astype(str)+ "_" + unique_df["other_pct_thresh"].astype(str)
#         unique_df["variant"] = unique_df.index
#         unique_df = unique_df.set_index("id")
#         all_unique_df.append(unique_df)
#
#     all_unique_df = pd.concat(all_unique_df)
#     all_unique_df["log2_n_cells"] = np.log2(all_unique_df["n_cells"]+1)
#     return all_unique_df
#
#
#
# def _objective_two_unique_vars_in_clone(all_unique_df, to_pivot=True):
#     if to_pivot:
#         if "id" in all_unique_df.columns:
#             df = all_unique_df.pivot(index="id", columns="variant", values="n_cells").fillna(0).astype(int)
#         else:
#             df = all_unique_df.reset_index().pivot(index="id", columns="variant", values="n_cells").fillna(0).astype(int)
#     else:
#         df = all_unique_df
#     vars_to_keep = df.loc[:,(df>0).sum()==1].columns # Variants with just 1 clone
#     clones_to_keep = df.loc[df.sum(axis=1)>1].index # Clones w >2 enriched variants
#     obj = 0
#     def cl_more_than_one(cl_ser):
#         curr = cl_ser[cl_ser > 0] # variants in a clone
#         # check if more than one unique variant for this clone
#         return sum([True if x in vars_to_keep else False for x in curr.index]) > 1
#
# #     def cl_more_than_one(cl_ser, vars_to_keep):
# #         curr = cl_ser[cl_ser > 0] # variants in a clone
# #         # check if more than one unique variant for this clone
# #         return sum([True if x in vars_to_keep else False for x in curr.index]) > 1
#     #obj = sum(df.loc[clones_to_keep].apply(lambda x: x(lambda y):, axis=0))
#     obj = sum(df.loc[clones_to_keep].apply(cl_more_than_one, axis=1))
#     return obj
#
#
# def _objectives(data):
#     all_unique_df = data["all_unique_df"]
#     #print('all_unique_df', all_unique_df.head)
#     ic('all_unique_df', all_unique_df.shape)
#     obj_max_nce_over_ncl = 0
#     obj_max_nce_over_nce = 0
#     obj_cl_over_ncl = 0
#     obj_nvars = 0
#     if len(all_unique_df) == 0:
#         #print('all 0', all_unique_df.columns)
#         return {x:(-1*np.inf) for x in objectives_l} # return score of 0 since all positive values
#     obj_d = all_unique_df.iloc[0]["pct_thresh"]
#     obj_e = all_unique_df.iloc[0]["other_pct_thresh"]
#     for v, v_df in all_unique_df.groupby("variant"):
#         ic(v)
#         max_ncells = max(v_df["n_cells"])
#         n_clones = len(set(v_df["clone"].values))
#         obj_max_nce_over_ncl += max_ncells/n_clones
#         obj_max_nce_over_nce += max_ncells/v_df["n_cells"].sum()
#
#         if n_clones != 0:
#             obj_cl_over_ncl += 1/n_clones
#             obj_nvars += 1
#
#     # calculate objective number of clones with more than one unique variant
#     obj_nclones_more_than_one_unique =  _objective_two_unique_vars_in_clone(all_unique_df, to_pivot=True)
#
#     objectives = {"variants_with_clone_norm_by_1_over_nclones_with_variant":obj_cl_over_ncl,
#                   "max_clone_ncells_over_nclones":obj_max_nce_over_ncl,
#                   "max_clone_ncells_over_ncells":obj_max_nce_over_nce,
#                   "pct_thresh":obj_d,"other_pct_thresh":obj_e,
#                    "n_vars":obj_nvars, "obj_nclones_more_than_one_unique": obj_nclones_more_than_one_unique}
#     return objectives
#
# def _constraints(solution):
#     #if solution["pct_thresh"] < solution["other_pct_thresh"]:
#     if "coverage_thresh" not in solution:
#         return None
#     if solution["af_thresh"]*solution["coverage_thresh"] >= 2:
#         return True
#     else:
#         return False
#
#
# def evaluate_series(individual_ser, AF_df, DP_df, curr_labels, return_data=False):
#     params = individual_ser.to_dict()
#     #print('params', params)
#     #solution = {"pct_thresh": individual[0], "af_thresh":individual[1],  "other_pct_thresh": individual[2]}
#     data = {"AF_df": AF_df, "DP_df":DP_df, "curr_labels":curr_labels}
#     all_unique_df = get_clones_unique_variants(params, data)
#     data["all_unique_df"] = all_unique_df
#     eval_out = _objectives(data)
#     if return_data:
#         return pd.Series(eval_out), data
#     else:
#         return pd.Series(eval_out)


# def set_multi(results, weights):
#     print(results.shape)
#     # first normalize results for each column to sum to 1
#     objs_total = results.replace([-np.inf, np.inf], np.nan).sum(axis=0)
#     print('objs_total', objs_total.head())
#     results_norm = results.apply(lambda x: x / objs_total.loc[x.name],
#                                  axis=0)
#
#     results_norm["multi"] = (weights * results_norm).sum(axis=1)
#     return results_norm.sort_values(by="multi")[::-1]
#
#
# def set_multi_rank(results, weights):
#     if "multi" in results.columns:  # in case multi was added before
#         rank_results = results.drop("multi", axis=1).rank(
#             na_option='top')
#     else:
#         rank_results = results.rank(na_option='top')
#     rank_results["multi"] = (weights * rank_results).sum(axis=1)
#     return rank_results.sort_values(by="multi")[::-1]

# def draw_heatmap(*args, **kwargs):
#     data = kwargs.pop('data')
#     # param_names = kwargs.pop('param_names', None)
#     print(data.shape)
#     d = data.pivot(index=args[1], columns=args[0],
#                    values=args[2]).fillna(0)
#
#     # get all clones and variants
#     d_full = pd.DataFrame(index=clones_order, columns=variants_order)
#     d_full.loc[:, :] = 0
#     d_full.loc[d.index, d.columns] = d
#     d_full = d_full.astype(float)
#
#     # get cluster results with jaccard
#     g = sns.clustermap(d_full, method="single")
#     inds = g.dendrogram_row.dendrogram["leaves"]
#     cols = g.dendrogram_col.dendrogram["leaves"]
#     plt.close(g.fig)
#
#     sns.heatmap(d_full.iloc[inds, cols],
#                 cbar_kws=dict(orientation="vertical"), **kwargs)
#     plt.title(data["params_multi"].values[0])
#     return
