import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def calc_cell_pairs(cell_ser, clones_series):
    curr_id = cell_ser.name
    print('curr_id', curr_id)
    print('clones_series', clones_series.head())
    curr_clone = clones_series.astype(object).fillna("").loc[curr_id]
    print('curr_clone', curr_clone)
    print('here')
    if curr_clone == "" or curr_clone is None : #np.isnan(curr_clone):
        cell_ser.loc[clones_series.loc[clones_series.isnull()].index] = -2
    else:
        cell_ser.loc[:] = (curr_clone == clones_series.loc[cell_ser.index])
        cell_ser.loc[clones_series.loc[clones_series.isnull()].index] = -1
    return cell_ser


def wrap_calc_cell_pairs(all_methods_df):
    donor_all_methods_pairs_df = {}
    for d, curr_methods_df in all_methods_df.groupby("donor"):
        print('donor', d)
        donor_all_methods_pairs_df[d] = {}
        # d = 0
        # curr_methods_df = all_methods_df.loc[all_methods_df["donor"]==d]
        # curr_donor_all_methods_pairs_df = {}

        for i in curr_methods_df.drop(["condition", "donor"],
                                      axis=1).columns:
            print('i', i)
            cell_pairs_df = pd.DataFrame(index=curr_methods_df.index,
                                         columns=curr_methods_df.index)
            donor_all_methods_pairs_df[d][i] = cell_pairs_df.apply(
                calc_cell_pairs, args=(curr_methods_df[i],), axis=1)

    return donor_all_methods_pairs_df


def calc_meth_overlap(meth_a, meth_b, use_na=False):
    n_T_T = ((meth_a == True) & (meth_b == True)).values.sum()
    n_F_F = ((meth_a == False) & (meth_b == False)).values.sum()
    n_F_T = ((meth_a == False) & (meth_b == True)).values.sum()
    n_T_F = ((meth_a == True) & (meth_b == False)).values.sum()

    n_nas_mone = ((meth_a == -1) & (meth_b == -1)).values.sum()
    n_nas_mtwo = ((meth_a == -2) & (meth_b == -2)).values.sum()

    n_nas_both = ((meth_a < 0) & (meth_b < 0)).values.sum()
    n_nas_mis = (((meth_a < 0) & (meth_b >= 0)) | (
                (meth_a >= 0) & (meth_b < 0))).values.sum()

    n_together = n_T_T + n_F_F + n_F_T + n_T_F

    # out_meth = meth_a == meth_b

    if not use_na:
        meth_a = meth_a.drop(meth_a < 0)
        meth_b = meth_b.drop(meth_b < 0)

        inds = set(meth_a.index).intersection(set(meth_b.index))
        meth_a = meth_a.loc[inds]
        meth_b = meth_b.loc[inds]
    # jaccard(meth_a, meth_b)
    n_T_T_norm = n_T_T / (n_T_T + n_T_F)  # n_together
    n_T_F_norm = n_T_F / (n_T_T + n_T_F)
    n_F_T_norm = n_F_T / (n_F_F + n_F_T)
    n_F_F_norm = n_F_F / (n_F_F + n_F_T) # n_together
    return {"n_T_T": n_T_T, "n_F_F": n_F_F, "n_F_T": n_F_T,
            "n_T_F": n_T_F,

            "n_T_T_norm": n_T_T_norm, "n_F_F_norm": n_F_F_norm,
            "n_F_T_norm": n_F_T_norm, "n_T_F_norm": n_T_F_norm,

            "n_nas_mone": n_nas_mone, "n_nas_mtwo": n_nas_mtwo,
            "n_nas_both": n_nas_both, "n_nas_mis": n_nas_mis,
            "n_together": n_together, "n_agree": n_T_T + n_F_F}


def get_title(fname, params_files):
    print('fname', fname)
    curr_file = fname.replace("/concat", "").replace("/cells_meta.tsv", "")
    #print('curr_file', curr_file)
    params = params_files.loc[curr_file].iloc[0]
    #print('params', params[["method", "nclonelist", "resolution"]])
    if params['method'] == 'knn':
        meth_params = params['resolution']
    if params['method'] == 'vireo':
        meth_params = params['nclonelist']
    return f"Variants method: {params['variants']}\n clones method: {params['method']}\nClones params: {meth_params}"

import itertools

def compare_methods(all_methods_pairs_df):
    meth_df = pd.DataFrame(
        columns=["donor", "m1", "m2", "n_agree", "n_T_T", "n_F_F",
                 "n_F_T", "n_T_F", "n_T_T_norm", "n_F_F_norm",
                 "n_F_T_norm", "n_T_F_norm", "n_nas_mone", "n_nas_mtwo",
                 "n_nas_both", "n_nas_mis", "n_together"])

    curr_d_methods = list(all_methods_pairs_df.keys())
    curr_d_pairs = list(itertools.product(curr_d_methods, repeat=2))
    for curr_pair in curr_d_pairs:
        if curr_pair[0] == curr_pair[1]:
            continue
        a = all_methods_pairs_df[curr_pair[0]]
        b = all_methods_pairs_df[curr_pair[1]]
        curr_out = calc_meth_overlap(a, b)
        #curr_out["donor"] = d
        curr_out["m1"] = curr_pair[0]
        curr_out["m2"] = curr_pair[1]
        meth_df = meth_df.append(pd.DataFrame(curr_out, index=[
            f"m1{curr_pair[0]}_m2{curr_pair[1]}"]))  # meth_df.loc[curr_pair]
    meth_df["n_agree_norm"] = meth_df["n_agree"] / meth_df["n_together"]
    return meth_df


sns.set(font_scale=1.2)
def draw_heatmap(*args, **kwargs):
    data = kwargs.pop('data')

    metric = args[0] #kwargs.pop(args[0])
    print('metric', metric)
    d = data.pivot(index="m1", columns="m2", values=[args[0]])
    #d = data.pivot(index=args[1], columns=args[0], values=args[2])
    sns.heatmap(d, **kwargs)
    return

def donors_heatmap(df, metric, title=""):
    df[metric] = df[metric].astype('float')
    g = sns.FacetGrid(df, col="donor", col_wrap=2, height=12, aspect=1)
    g.map_dataframe(draw_heatmap, metric, cbar=True, square=True)
    g.figure.suptitle(title)#"Methods comparison: #pairs in shared clusters clusters")
    plt.tight_layout()
    return


def donors_agg_heatmap(df, metric, title="", to_log=False):
    df = df.copy()
    if to_log:
        df[metric] = np.log10(df[metric]+1)
    df_agg = df.groupby(['m1','m2']).mean().reset_index()
    f = plt.figure(figsize=(12,12), dpi=300)
    draw_heatmap(metric, data=df_agg)
    plt.title(title)#"Methods comparison: #pairs in shared clusters clusters")
    plt.tight_layout()
    return