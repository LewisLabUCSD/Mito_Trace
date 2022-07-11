""" Snakemake Script. Don't need to import snakemake."""

import pandas as pd
from os.path import dirname, join
import numpy as np
#import click
#import snakemake
from src.utils.data_io import af_to_vireo
from icecream import ic
import seaborn as sns
import matplotlib.pyplot as plt


def run_impute(af_df, dp_df, cov_thresh, cell_pct_cov_thresh):
    dp_low_cov = dp_df < cov_thresh
    vars_pct_cov = dp_low_cov.sum()/dp_df.shape[0]
    filt_vars_pct_cov = vars_pct_cov[(1-vars_pct_cov)>cell_pct_cov_thresh]
    dp_low_cov = dp_low_cov.reset_index().melt(id_vars=["Cell"], var_name="Variant",value_name="low")
    print("before filt", dp_low_cov.shape)
    dp_low_cov = dp_low_cov[dp_low_cov["low"]]
    print("after filt, low cov", dp_low_cov.shape)
    dp_var_group = dp_low_cov.groupby("Variant")

    imp_vars = set(dp_low_cov["Variant"].values)
    imp_df = af_df.copy()
    inds_imp = {}
    median_dict = {}
    for v in imp_vars:
        if v in filt_vars_pct_cov.index:
            curr_lowc = dp_var_group.get_group(v)["Cell"].values
            #print(sum(af_df.index.isin(curr_lowc)))
            median_v = af_df.loc[~(af_df.index.isin(curr_lowc)), v].mean()
            # print('variant', v)
            # print('median', median_v )
            #print(sum(imp_df.index.isin(curr_lowc)))
            imp_df.loc[imp_df.index.isin(curr_lowc), v] = median_v
            inds_imp[v] = curr_lowc
            median_dict[v] = median_v
        # else:
        #     print('variant low coverage', v)
    return imp_df, inds_imp, median_dict


def plot_impute(imp_df, af_df, dp_df, cov_thresh, cell_pct_cov_thresh, median_d, f_save=None, include_orig=False):
    #g_pile = sns.clustermap(imp_df)

    #imp_df
    col_meta = pd.DataFrame(median_d, index=["median"]).transpose()

    dp_df = dp_df.loc[af_df.index, af_df.columns]
    for v in imp_df.columns:
        if v not in col_meta.index:
            col_meta.loc[v, "median"] = 0 #np.nan


    from mplh import cluster_help as ch
    g_pile=ch.plot_cluster(df=imp_df,
                    col_meta=col_meta,
                    col_clr_schemes="sequential")
    plt.suptitle(f"imputation with donor cov_thresh {cov_thresh} cell_pct_cov_thresh {cell_pct_cov_thresh}")
    inds = g_pile.dendrogram_row.dendrogram['leaves']
    cols = g_pile.dendrogram_col.dendrogram["leaves"]

    g_dp = ch.plot_cluster(df=np.log2(1+dp_df.iloc[inds, cols]),
                    to_row_clust=False, to_col_clust=False,
                    col_meta=col_meta, col_clr_schemes="sequential")
    #sns.clustermap(np.log2(1+dp_df.iloc[inds, cols]), row_cluster=False, col_cluster=False)
    plt.suptitle(f"no imputation cov_thresh {cov_thresh} cell_pct_cov_thresh {cell_pct_cov_thresh}")
    if f_save is not None:
        g_pile.fig.savefig(f_save+".af.tsv")
        g_dp.fig.savefig(f_save+".dp.tsv")

    if include_orig:
        sns.clustermap(af_df.iloc[inds, cols], row_cluster=False, col_cluster=False)
        plt.suptitle(f"no imputation cov_thresh {cov_thresh} cell_pct_cov_thresh {cell_pct_cov_thresh}")

    return g_pile


def impute_and_plot(af_df, dp_df, cov_thresh, cell_pct_cov_thresh, f_save=None, include_orig=False):
    imp_df, inds_imp, median_dict = run_impute(af_df, dp_df, cov_thresh, cell_pct_cov_thresh)
    g_pile = plot_impute(imp_df, af_df, dp_df, cov_thresh, cell_pct_cov_thresh, median_dict, f_save=f_save, include_orig=include_orig)
    print(f"n vars imputed {len(median_dict)}")
    return {"impute_af": imp_df, "figure": g_pile, "median_d": median_dict,
     "imputed_indices": inds_imp}

