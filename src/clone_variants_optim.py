import pandas as pd
import numpy as np
from icecream import ic
import seaborn as sns
import matplotlib.pyplot as plt

ic.disable()


def get_unique_variants(cln_af, other_af, pct_thresh, af_thresh,
                        other_pct_thresh):
    """ gets the distinct variants in a clone.
    """
    n_thresh = pct_thresh * cln_af.shape[1]
    n_oth_thresh = other_pct_thresh * other_af.shape[1]
    bin_cln = cln_af > af_thresh
    bin_other = other_af > af_thresh
    cells_above = bin_cln.sum(axis=1)
    pct_cells_above = cells_above / bin_cln.shape[1]
    up_vars = bin_cln.loc[cells_above > n_thresh].index
    cells_other_above = bin_other.sum(axis=1)
    pct_cells_other_above = cells_other_above / bin_other.shape[1]
    up_oth_vars = bin_other.loc[cells_other_above > n_oth_thresh].index
    uniq_vars = list(set(up_vars) - set(up_oth_vars))
    out = pd.DataFrame(index=uniq_vars, data={
        "n_cells": cells_above.loc[uniq_vars].values,
        "n_other_cells": cells_other_above.loc[uniq_vars].values,
        "pct_above": pct_cells_above,
        "pct_other_above": pct_cells_other_above})
    out["pct_thresh"] = pct_thresh
    out["af_thresh"] = af_thresh
    out["other_pct_thresh"] = other_pct_thresh
    return out


def get_clones_unique_variants(solution, data):
    all_unique_df = []
    pct_thresh, af_thresh, other_pct_thresh = solution["pct_thresh"], \
                                              solution["af_thresh"], \
                                              solution[
                                                  "other_pct_thresh"]  # solution[0], solution[1], solution[2]
    curr_labels = data["curr_labels"]
    AF_df = data["AF_df"]
    DP_df = data["DP_df"]
    for cln, val in curr_labels.groupby("name"):
        ic(cln)
        cln_af = AF_df.loc[:, val.index]
        other_af = AF_df.loc[:, curr_labels.drop(val.index).index]
        curr_dp = DP_df.loc[:, val.index]
        curr_labs = curr_labels[curr_labels.index.isin(cln_af.columns)]
        ic(cln_af.shape)
        unique_df = get_unique_variants(cln_af, other_af, pct_thresh,
                                        af_thresh, other_pct_thresh)
        unique_df["clone"] = cln
        unique_df["id"] = unique_df["clone"] + "_" + unique_df[
            "pct_thresh"].astype(str) + "_" + unique_df[
                              "af_thresh"].astype(str) + "_" + \
                          unique_df["other_pct_thresh"].astype(str)
        unique_df["variant"] = unique_df.index
        unique_df = unique_df.set_index("id")
        all_unique_df.append(unique_df)

    all_unique_df = pd.concat(all_unique_df)
    all_unique_df["log2_n_cells"] = np.log2(
        all_unique_df["n_cells"] + 1)
    return all_unique_df


def _objective_two_unique_vars_in_clone(all_unique_df, to_pivot=True):
    if to_pivot:
        if "id" in all_unique_df.columns:
            df = all_unique_df.pivot(index="id", columns="variant",
                                     values="n_cells").fillna(0).astype(
                int)
        else:
            df = all_unique_df.reset_index().pivot(index="id",
                                                   columns="variant",
                                                   values="n_cells").fillna(
                0).astype(int)
    else:
        df = all_unique_df
    vars_to_keep = df.loc[:, (
                                         df > 0).sum() == 1].columns  # Variants with just 1 clone
    clones_to_keep = df.loc[
        df.sum(axis=1) > 1].index  # Clones w >2 enriched variants
    obj = 0

    def cl_more_than_one(cl_ser):
        curr = cl_ser[cl_ser > 0]  # variants in a clone
        # check if more than one unique variant for this clone
        return sum([True if x in vars_to_keep else False for x in
                    curr.index]) > 1

    #     def cl_more_than_one(cl_ser, vars_to_keep):
    #         curr = cl_ser[cl_ser > 0] # variants in a clone
    #         # check if more than one unique variant for this clone
    #         return sum([True if x in vars_to_keep else False for x in curr.index]) > 1
    # obj = sum(df.loc[clones_to_keep].apply(lambda x: x(lambda y):, axis=0))
    obj = sum(df.loc[clones_to_keep].apply(cl_more_than_one, axis=1))
    return obj


def _objectives(data, objectives_l):
    all_unique_df = data["all_unique_df"]
    # print('all_unique_df', all_unique_df.head)
    ic('all_unique_df', all_unique_df.shape)
    obj_max_nce_over_ncl = 0
    obj_max_nce_over_nce = 0
    obj_cl_over_ncl = 0
    obj_nvars = 0
    if len(all_unique_df) == 0:
        # print('all 0', all_unique_df.columns)
        return {x: (-1 * np.inf) for x in
                objectives_l}  # return score of 0 since all positive values
    obj_d = all_unique_df.iloc[0]["pct_thresh"]
    obj_e = all_unique_df.iloc[0]["other_pct_thresh"]
    for v, v_df in all_unique_df.groupby("variant"):
        ic(v)
        max_ncells = max(v_df["n_cells"])
        n_clones = len(set(v_df["clone"].values))
        obj_max_nce_over_ncl += max_ncells / n_clones
        obj_max_nce_over_nce += max_ncells / v_df["n_cells"].sum()

        if n_clones != 0:
            obj_cl_over_ncl += 1 / n_clones
            obj_nvars += 1

    # calculate objective number of clones with more than one unique variant
    obj_nclones_more_than_one_unique = _objective_two_unique_vars_in_clone(
        all_unique_df, to_pivot=True)

    objectives = {
        "variants_with_clone_norm_by_1_over_nclones_with_variant": obj_cl_over_ncl,
        "max_clone_ncells_over_nclones": obj_max_nce_over_ncl,
        "max_clone_ncells_over_ncells": obj_max_nce_over_nce,
        "pct_thresh": obj_d, "other_pct_thresh": obj_e,
        "n_vars": obj_nvars,
        "obj_nclones_more_than_one_unique": obj_nclones_more_than_one_unique}
    return objectives


def _constraints(solution):
    # if solution["pct_thresh"] < solution["other_pct_thresh"]:
    if "coverage_thresh" not in solution:
        return None
    if solution["af_thresh"] * solution["coverage_thresh"] >= 2:
        return True
    else:
        return False


def evaluate_series(individual_ser, AF_df, DP_df, curr_labels,
                    return_data=False,
                    objectives_l=("variants_with_clone_norm_by_1_over_nclones_with_variant",
                                  "max_clone_ncells_over_nclones",
                                  "max_clone_ncells_over_ncells",
                                  "pct_thresh","other_pct_thresh",
                                  "n_vars",
                                  "obj_nclones_more_than_one_unique")):
    params = individual_ser.to_dict()
    # print('params', params)
    # solution = {"pct_thresh": individual[0], "af_thresh":individual[1],  "other_pct_thresh": individual[2]}
    data = {"AF_df": AF_df, "DP_df": DP_df, "curr_labels": curr_labels}
    all_unique_df = get_clones_unique_variants(params, data)
    data["all_unique_df"] = all_unique_df
    eval_out = _objectives(data, objectives_l=objectives_l)
    if return_data:
        return pd.Series(eval_out), data
    else:
        return pd.Series(eval_out)


def set_multi_rank(results, weights):
    if "multi" in results.columns:  # in case multi was added before
        rank_results = results.drop("multi", axis=1).rank(
            na_option='top')
    else:
        rank_results = results.rank(na_option='top')
    rank_results["multi"] = (weights * rank_results).sum(axis=1)
    return rank_results.sort_values(by="multi")[::-1]


def set_multi(results, weights):
    print(results.shape)
    # first normalize results for each column to sum to 1
    objs_total = results.replace([-np.inf, np.inf], np.nan).sum(axis=0)
    print('objs_total', objs_total.head())
    results_norm = results.apply(lambda x: x / objs_total.loc[x.name],
                                 axis=0)

    results_norm["multi"] = (weights * results_norm).sum(axis=1)
    return results_norm.sort_values(by="multi")[::-1]

#
#
#
# heatmap_input = all_df[["n_cells", "variant"]].reset_index().pivot(
#     index="id", columns="variant", values="n_cells").fillna(0).astype(
#     int)
# meta_df = all_df[
#     ["af_thresh", "other_pct_thresh", "pct_thresh", "clone"]]
# meta_df = meta_df.loc[~(meta_df.index.duplicated())]
# meta_df = meta_df.sort_values(
#     ["af_thresh", "pct_thresh", "other_pct_thresh", "clone"])
# heatmap_input = heatmap_input.loc[meta_df.index]
#
# # Get the variants based on total number of cells across parameters
# heatmap_input = heatmap_input.loc[:,
#                 heatmap_input.sum().sort_values()[::-1].index]
# variants_order = heatmap_input.columns
# clone_sums = meta_df.groupby("clone").apply(clone_sum)
# clone_sums = clone_sums.loc[:,
#              clone_sums.sum().sort_values()[::-1].index]
# clones_order = clone_sums.index
#

# Get the clones based on total number of cells across parameters
def clone_sum(val, heatmap_input):
    # print(val)
    return (heatmap_input.loc[val.index].sum())


def prep_long_heatmap(all_df):
    heatmap_input = all_df[["n_cells", "variant"]].reset_index().pivot(
        index="id", columns="variant", values="n_cells").fillna(
        0).astype(int)
    meta_df = all_df[
        ["af_thresh", "other_pct_thresh", "pct_thresh", "clone"]]
    meta_df = meta_df.loc[~(meta_df.index.duplicated())]
    meta_df = meta_df.sort_values(
        ["af_thresh", "pct_thresh", "other_pct_thresh", "clone"])
    heatmap_input = heatmap_input.loc[meta_df.index]

    # Get the variants based on total number of cells across parameters
    heatmap_input = heatmap_input.loc[:,
                    heatmap_input.sum().sort_values()[::-1].index]
    variants_order = heatmap_input.columns
    clone_sums = meta_df.groupby("clone").apply(clone_sum,
                                                heatmap_input)
    clone_sums = clone_sums.loc[:,
                 clone_sums.sum().sort_values()[::-1].index]
    clones_order = clone_sums.index
    return clones_order, variants_order, heatmap_input

def params_to_str(ser, param_names):
    name = ""
    for p in param_names:
        name = name + f"{p}={ser[p]:.3f}\n"
    return name

def params_and_multi_str(ser):
    param_str = ser["params"]
    name = f"params:\n{param_str.strip()}\nObjective score={ser['multi_obj']} (want to maximize)"
    return name


def draw_heatmap(*args, **kwargs):
    data = kwargs.pop('data')
    clones_order = kwargs.pop('clones_order')
    variants_order = kwargs.pop('variants_order' )

    #title_col = kwargs.pop("title_col", "params_multi")
    #print('title_col', title_col)
    # param_names = kwargs.pop('param_names', None)
    # print(data.shape)
    # print(data.head())
    d = data.pivot(index=args[1], columns=args[0],
                   values=args[2]).fillna(0)

    # get all clones and variants
    d_full = pd.DataFrame(index=clones_order, columns=variants_order)
    d_full.loc[:, :] = 0
    d_full.loc[d.index, d.columns] = d
    d_full = d_full.astype(float)

    # get cluster results
    g = sns.clustermap(d_full, method="single")
    inds = g.dendrogram_row.dendrogram["leaves"]
    cols = g.dendrogram_col.dendrogram["leaves"]
    plt.close(g.fig)

    sns.heatmap(d_full.iloc[inds, cols],xticklabels=True, yticklabels=True,
                cbar_kws=dict(orientation="vertical"), **kwargs)

    #print(data[title_col].values[0])
    #plt.gca().set_title(data[title_col].values[0])
    return
