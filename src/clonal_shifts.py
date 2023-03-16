import pandas as pd
import numpy as np
from os.path import join
from tqdm.notebook import tqdm

from scipy.stats import hypergeom, fisher_exact
from statsmodels.stats import multitest

import seaborn as sns
import matplotlib.pyplot as plt
from icecream import ic
import pickle

print('here')
def merge_hypergeom(a_enrich_df, b_enrich_df, a_name, b_name,
                    p_thresh=0.1, f_save=None,
                    title="Number of clones that are significant"):
    """ Combine two hypergeometric tests that have the same dimensions.
    Looks for samples (clones) that have significance both tests

    :param a_enrich_df:
    :param b_enrich_df:
    :param a_name:
    :param b_name:
    :param p_thresh:
    :param f_save:
    :param title:
    :return:
    """
    a_enrich_df.index.name = "name"
    b_enrich_df.index.name = "name"

    a_sig = a_enrich_df < p_thresh
    if a_sig.shape[0] == 0:
        a_sig = pd.DataFrame(columns=["name", a_name, "sig"])
    else:
        a_sig = a_sig.reset_index().melt(id_vars="name", value_name="sig",
                                         var_name=a_name)
    print('a_sig', a_sig.head())
    b_sig = b_enrich_df < p_thresh
    print(b_sig.shape)
    if b_sig.shape[0] == 0:
        b_sig = pd.DataFrame(columns=["name", b_name, "sig"])
    else:
        b_sig = b_sig.reset_index().melt(id_vars="name", value_name="sig",
                                         var_name=b_name)

    print('b_sig', b_sig.head())

    merged_df = pd.DataFrame(index=a_sig[a_name].unique(),
                             columns=b_sig[b_name].unique()).astype(object)
    merged_df.loc[:, :] = ""
    merged_count_df = pd.DataFrame(index=a_sig[a_name].unique(),
                                   columns=b_sig[
                                       b_name].unique()).astype("Int64")
    merged_count_df = merged_count_df.fillna(0)
    #print(b_sig.loc[b_sig["sig"]])
    b_sig = b_sig.loc[b_sig["sig"]].drop("sig", axis=1)
    a_sig = a_sig.loc[a_sig["sig"]].drop("sig", axis=1)

    for inp_ind, val in a_sig.iterrows():
        curr_clone = val["name"]
        curr_cult_df = b_sig.loc[b_sig["name"] == curr_clone]
        for noinp_ind, val2 in curr_cult_df.iterrows():
            merged_df.loc[val[a_name], val2[b_name]] = merged_df.loc[
                                                           val[a_name],
                                                           val2[
                                                               b_name]] + curr_clone + ";"
            merged_count_df.loc[val[a_name], val2[b_name]] += 1
    print('merge_df', merged_df)
    print('count', merged_count_df)
    merged_df = merged_df.apply(
        lambda ser: ser.apply(lambda x: x.strip(";")), axis=0)
    merged_df.index.name = a_name
    merged_df.columns.name = b_name

    f, ax = plt.subplots(figsize=(12, 12), dpi=300)
    try:
        sns.heatmap(merged_count_df.astype(int),
                    annot=merged_df, fmt="s")
        plt.ylabel(a_name)
        plt.xlabel(b_name)
        plt.title(title)
    except ValueError:  # raised if `y` is empty.
        print('No clones significant')
        plt.title(f"{title} no significant clones")
        pass

    if f_save is not None:
        plt.savefig(f_save + ".pdf")
        merged_df.to_csv(f_save + ".csv")
    return merged_df, merged_count_df


def get_groups(groups, clones=None, atac_cl=None, clone_col=None,
               atac_col=None):
    """ Gets the unique values for the columns based on the groups df

    :param groups: DataFrame with the atac_col and clone_col
    :param clones: if not None, skips getting the unique sample columns
    :param atac_cl: if not None, skips getting the unique feature columns
    :param clone_col: sample columns name
    :param atac_col: feature column name .
    :return:
    """
    if clones is None:
        if clone_col is None:
            raise ValueError("clone_col or clones cant be None")
        clones = groups[clone_col].unique()
    if atac_cl is None:
        if atac_col is None:
            raise ValueError("atac_col or atac_cl cant be None")
        atac_cl = groups[atac_col].unique()
    return clones, atac_cl


def run_hypergeo(groups, clones, atac_cl, atac_col, clone_col):
    """ Computes hypergeometric for each feature for each sample.

    :param groups:
    :param clones:
    :param atac_cl:
    :param atac_col:
    :param clone_col:
    :return:
    """
    # p(k,M,n,N) = (n choose k)((M-n)choose(N-k))/(MchooseN)
    # pmf(k, M, n, N) = choose(n, k) * choose(M - n, N - k) / choose(M, N),
    # for max(0, N - (M-n)) <= k <= min(n, N)

    # M: Total number of cells
    # n: Number of cells in the atac cluster (group population)
    # N: Number of cells in clone (the draw)
    # x: Number of cells in specific clone and cluster
    enrichment_df = pd.DataFrame(index=clones, columns=atac_cl,
                                 dtype=np.float128)

    M = groups["count"].sum()
    for cl in clones:
        for atac in atac_cl:
            # n = groups[groups[atac_col] == atac]["count"].sum()
            # N = groups[groups[clone_col] == cl]["count"].sum()
            N = groups[groups[atac_col] == atac]["count"].sum()
            n = groups[groups[clone_col] == cl]["count"].sum()

            x = groups[((groups[clone_col] == cl) & (
                        groups[atac_col] == atac))]["count"].sum()

            # rv = hypergeom(M, n, N)
            prb = 1 - hypergeom.cdf(x, M, n, N)
            enrichment_df.loc[cl, atac] = prb
    return enrichment_df


def pval_correct(enrichment_df, p_thresh):
    """ Does BH-FDR test for hypergeometric tests

    :param enrichment_df:
    :param p_thresh:
    :return:
    """
    reject, pvals_corrected, _, _ = multitest.multipletests(
        enrichment_df.values.flatten(), alpha=p_thresh, method="fdr_bh")

    nrows, ncols = enrichment_df.shape
    #print(enrichment_df.values.flatten())
    reject, pvals_corrected, _, _ = multitest.multipletests(
        enrichment_df.values.flatten(), alpha=p_thresh, method="fdr_bh")
    pvals_corrected = np.reshape(pvals_corrected, [nrows, ncols])
    return pvals_corrected


def create_enrichment(groups, atac_col, clone_col, p_thresh, clones=None, atac_cl=None, to_correct=True):
    """ Runs the hypergeometric test and does p-val correction if to_correct=True

    :param groups:
    :param atac_col:
    :param clone_col:
    :param p_thresh:
    :param clones:
    :param atac_cl:
    :param to_correct:
    :return:
    """
    clones, atac_cl = get_groups(groups, clones, atac_cl, clone_col,
                                 atac_col)
    enrichment_df = run_hypergeo(groups, clones, atac_cl, atac_col, clone_col)
    if to_correct:
        pvals_corrected = pval_correct(enrichment_df, p_thresh)
        bh_enrichment_df = enrichment_df.copy()
        bh_enrichment_df.loc[:, :] = pvals_corrected
        # print('bh_enrichment_df', bh_enrichment_df.head())
        return bh_enrichment_df
    return enrichment_df


def set_minimum_p(df, p_thresh):
    # sets the minimum p-value instead of 0
    min_p = (set((df.values).flatten()) - {0})
    if len(min_p) == 0:
        min_p_val = p_thresh-(p_thresh/10**4)
    else:
        min_p_val = min(min(min_p), p_thresh-(p_thresh/10**4))
    df[df == 0] = min_p_val  # Set to the next min, or p_thresh, whichever is smaller
    return df

def create_output_df(bh_enrichment_df, sizes, p_thresh, to_filt=True):
    """ Calculates significance with p_thresh and filters if to_filt

    :param bh_enrichment_df:
    :param sizes:
    :param p_thresh:
    :param to_filt:
    :return:
    """
    output_df = pd.DataFrame(index=sizes.index)
    # output_df = pd.DataFrame(index=bh_enrichment_df.index)
    output_df["significant clusters"] = ""
    output_df["size"] = sizes
    output_df["min_significance"] = None

    # for ind, val in bh_enrichment_df.loc[sizes.index].iterrows():
    for ind, val in bh_enrichment_df.iterrows():
        passed = val[val < p_thresh].index.values
        if len(passed) > 0:
            output_df.loc[ind, "significant clusters"] = ";".join(
                [str(x) for x in passed])
            output_df.loc[ind, "min_significance"] = min(
                val)  # sig_results.append((ind, passed))
    output_df.loc[:, bh_enrichment_df.columns] = bh_enrichment_df.loc[
        output_df.index]


    output_df = output_df.loc[~(output_df["min_significance"].isnull())]
    output_df = output_df.sort_values("min_significance")
    #output_df = output_df.sort_values("size", ascending=True)
    if to_filt:
        bh_enrichment_df[bh_enrichment_df > p_thresh] = 1

    bh_enrichment_df = set_minimum_p(bh_enrichment_df, p_thresh)
    return output_df, bh_enrichment_df


def pipeline_groups_hypergeo(groups, clones, atac_cl, sizes, p_thresh, atac_col, clone_col, to_filt=True, to_correct=True):
    """ Runs bh-corrected hypergeo for each clone-cluster combination

    :param groups:
    :param clones:
    :param atac_cl:
    :param sizes:
    :param p_thresh:
    :param atac_col:
    :param clone_col:
    :param to_filt:
    :param to_correct:
    :return:
    """
    bh_enrichment_df =  create_enrichment(groups, atac_col, clone_col,
                                          p_thresh, clones, atac_cl, to_correct=to_correct)
    output_df, bh_enrichment_df = create_output_df(bh_enrichment_df,
                                                   sizes=sizes, p_thresh=p_thresh, to_filt=to_filt)
    return output_df, bh_enrichment_df


def run(groups, sizes, atac_col, clone_col, p_thresh,
        clones=None, atac_cl=None, out_f=None, to_plot=True):
    if clones is None:
        clones = groups[clone_col].unique()
    if atac_cl is None:
        atac_cl = groups[atac_col].unique()
    clones, atac_cl = get_groups(groups, clones, atac_cl, clone_col,atac_col)
    output_df, bh_enrichment_df = pipeline_groups_hypergeo(groups,
                                                           clones,
                                                           atac_cl,
                                                           sizes,
                                                           atac_col=atac_col,
                                                           clone_col=clone_col,
                                                           p_thresh=p_thresh,
                                                           )
    if to_plot:
        if output_df.shape[0] == 0:
            f = plt.figure()
            plt.heatmap(-np.log10(bh_enrichment_df.fillna(1)),row_cluster=False)
            plt.title("-log10 p-value (No groups were significant)")
        else:
            g = sns.clustermap(
                -np.log10(bh_enrichment_df.loc[output_df.index].fillna(1)))
            g.ax_heatmap.set(xlabel="Cluster ID")
            g.ax_cbar.set(title="-log10 p-value")

        if out_f is not None:
            plt.savefig(out_f,dpi=300, bbox_inches = "tight")
        # plot just the counts
        groups["log2_count"] = np.log2(groups["count"] + 1)
        f = plt.figure()
        sns.clustermap(
            groups.pivot(index="cluster_labels", columns="den_clust",
                              values="log2_count").fillna(0))
        if out_f is not None:
            plt.savefig(out_f + "groups_counts.png",dpi=300, bbox_inches = "tight")
            plt.savefig(out_f + "groups_counts.pdf",dpi=300, bbox_inches = "tight")
    return output_df, bh_enrichment_df


def wrap_create_enrichment(ar, groups, atac_col, clone_col, p_thresh, clones, atac_cl, to_p_correct):
    out = []
    cell_pairs = []

    for ind, val in groups.iterrows():
        # if val["count"] > 0
        cell_pairs.extend(
            [[val[atac_col], val[clone_col]]] * int(val["count"]))
    cell_pairs = np.array(cell_pairs)

    for i in ar:
        #ic(i)
        #curr_sim = groups.copy()
        #curr_sim[clone_col] = curr_sim[clone_col].sample(n=groups.shape[0]).values
        cell_pairs[:, 1] = np.random.choice(cell_pairs[:, 1],
                                            size=cell_pairs.shape[0],
                                            replace=False)
        curr_sim= pd.DataFrame(cell_pairs, columns=[atac_col,clone_col]).groupby([atac_col, clone_col]).size().reset_index().rename({0: "count"}, axis=1)

        #print('curr_sim', curr_sim["count"].sum())
        out.append(create_enrichment(curr_sim, atac_col, clone_col, p_thresh,
                          clones=clones, atac_cl=atac_cl, to_correct=to_p_correct))
        # print('after')
        # print(curr_sim.head())
    return out


def shuffle_hypergeo(groups, atac_col, clone_col, p_thresh, clones, atac_cl, to_p_correct=True, n_shuffle=1000,
                     to_parallel=True, n_cpus=4):
    # expand groups to shuffle - convert to long vector
    # clones, atac_cl = get_groups(groups, clones, atac_cl, clone_col,
    #            atac_col)
    # all_out = pd.DataFrame(index=clones, columns=atac_cl)

    if to_parallel:
        #from numpanpar import parallel_ar as parar
        from src.utils.parallel_helper import parallel_ar as parar
        all_out = parar(np.arange(n_shuffle), func=wrap_create_enrichment,
                        func_args=(groups, atac_col, clone_col, p_thresh, clones, atac_cl, to_p_correct),num_processes=n_cpus, to_flat=False)
    else:
        all_out = wrap_create_enrichment(np.arange(n_shuffle), groups, atac_col, clone_col, p_thresh, clones, atac_cl, to_p_correct)
    return all_out


def get_shuffle_results(shuffle, clones, clone_map):
    # a. Get the min for each run
    # b. Get the min for each clone
    # c. Get all values across all runs
    # d. Get all values for each clone separatelyt

    global_min = [i.min().min() for i in shuffle]
    #ic(len(global_min))
    clone_min = {}
    for curr_ind in clones:
        clone_min[curr_ind] = [i[clone_map[curr_ind]].min() for i in
                               shuffle]
    #ic(len(clone_min))

    clone_min = {}
    for curr_ind in clones:
        clone_min[curr_ind] = [i[clone_map[curr_ind]].min() for i in
                               shuffle]
    #ic(len(clone_min))

    global_all = []
    for i in shuffle:
        global_all.extend(i.flatten())
    len(global_all)

    clone_all = {}
    for curr_ind in clones:
        clone_all[curr_ind] = []
        for i in shuffle:
            clone_all[curr_ind].extend(i[clone_map[curr_ind]])
        #ic(len(clone_all[curr_ind]))
    #ic(len(clone_all))

    return global_min, clone_min, global_all, clone_all


def plot_glob_all(bh_enrichment_df, global_all, p_thresh, atac_col, out_f=None):
    shuffle_pval_df = bh_enrichment_df.copy().applymap(lambda x: ((x>np.array(global_all)).sum())/len(global_all))
    #shuffle_pval_df = ((bh_enrichment_df) < global_all).sum() / len(global_all)
    curr_thresh = np.percentile(global_all, 100 * p_thresh)
    f, ax = plt.subplots()
    sns.histplot(global_all, ax=ax, stat='probability', binwidth=0.05, binrange=(0,1))
    sns.histplot(bh_enrichment_df.values.flatten(), color="green", stat='probability', alpha=0.4, ax=ax,  binwidth=0.05, binrange=(0,1))
    plt.axvline(curr_thresh, color='r', alpha=0.4)
    group_vs_shuffle = bh_enrichment_df[bh_enrichment_df < curr_thresh]
    sig_pairs = group_vs_shuffle.stack().reset_index().rename({"level_1": atac_col, 0: "p_val"}, axis=1)
    print(f"Number of groups below p-val significance: {(sig_pairs.shape[0])}")
    plt.title(f"p-value across all simulations\ncutoff: {p_thresh}; nsig: {(sig_pairs.shape[0])}")
    sig_pairs = sig_pairs.sort_values("p_val")
    if out_f is not None:
        plt.savefig(out_f,dpi=300, bbox_inches = "tight")
    return sig_pairs, shuffle_pval_df, curr_thresh


def plot_glob_min(bh_enrichment_df, global_min, p_thresh, atac_col, out_f=None):
    shuffle_pval_df = bh_enrichment_df.copy().applymap(
        lambda x: ((x > np.array(global_min)).sum()) / len(global_min))
    #shuffle_pval_df = ((bh_enrichment_df) < global_min).sum() / len(global_min)
    curr_thresh = np.percentile(global_min, 100 * p_thresh)
    f, ax = plt.subplots()
    sns.histplot(global_min, ax=ax, stat='probability',  binwidth=0.05, binrange=(0,1))
    sns.histplot(bh_enrichment_df.values.flatten(), color="green",  binwidth=0.05, binrange=(0,1), stat='probability', alpha=0.4, ax=ax)
    plt.axvline(curr_thresh, color='r')
    group_vs_shuffle = bh_enrichment_df[bh_enrichment_df < curr_thresh]
    sig_pairs = group_vs_shuffle.stack().reset_index().rename(
        {"level_1": atac_col, 0: "p_val"}, axis=1)
    print(f"Number of groups below p-val significance: {(sig_pairs.shape[0])}")
    sig_pairs = sig_pairs.sort_values("p_val")
    plt.title(f"Minimum p-value across all simulations\ncutoff: {p_thresh}; nsig: {(sig_pairs.shape[0])}")
    if out_f is not None:
        plt.savefig(out_f,dpi=300, bbox_inches = "tight")
    return sig_pairs, shuffle_pval_df, curr_thresh


## clone_min
def plot_clone_min(bh_enrichment_df, clone_min, p_thresh, out_f=None):
    clone_min_sig = {}

    g = sns.FacetGrid(pd.DataFrame(clone_min).melt(), col="variable",
                      col_wrap=4, sharex=False, sharey=False)
    g.map_dataframe(sns.histplot, x="value", stat='probability', binwidth=0.05, binrange=(0,1),
                    color='blue', alpha=0.6)

    shuffle_pval_df = bh_enrichment_df.copy()
    for i, ind in enumerate(bh_enrichment_df.index):
        shuffle_pval_df.loc[ind] = bh_enrichment_df.loc[ind].apply(
            lambda x: (x > np.array(clone_min[ind])).sum() / len(
                clone_min[ind]))
        #shuffle_pval_df.loc[ind] = ((bh_enrichment_df.loc[ind])< clone_min[ind]).sum()/len(clone_min[ind])
        curr_thresh = np.percentile(clone_min[ind], 100 * p_thresh)
        group_vs_shuffle = bh_enrichment_df.loc[ind] < curr_thresh
        sig_pairs = group_vs_shuffle[group_vs_shuffle]
        clone_min_sig[ind] = sig_pairs
        g.axes_dict[ind].axvline(curr_thresh, color='red', linewidth=4)
        sns.histplot(bh_enrichment_df.loc[ind], stat='probability', binwidth=0.05, binrange=(0,1),
                     color='green', alpha=0.4, ax=g.axes_dict[ind])
    if out_f is not None:
        plt.savefig(out_f,dpi=300, bbox_inches = "tight")
    return clone_min_sig, shuffle_pval_df, curr_thresh


def plot_clone_all(bh_enrichment_df, clone_all, p_thresh, out_f=None):
    f, axs = plt.subplots(
        nrows=int(np.ceil(bh_enrichment_df.shape[0] / 4)), ncols=4,
        figsize=(12, 12), sharex=False, sharey=False, squeeze=False)
    print('axs', len(axs))
    print('bh shape', bh_enrichment_df.shape)
    clone_min_sig = {}
    shuffle_pval_df = bh_enrichment_df.copy()
    for i, ind in enumerate(bh_enrichment_df.index):
        #print('i,ind', i, ind)
        curr_row, curr_col = int(np.floor(i / 4)), i % 4
        #print('curr row', curr_row, 'curr_col', curr_col)
        data = clone_all[ind]
        #print('data',data)

        shuffle_pval_df.loc[ind] = bh_enrichment_df.loc[ind].apply(lambda x: (x > np.array(clone_all[ind])).sum()/len(clone_all[ind]))

        curr_thresh = np.percentile(clone_all[ind], 100 * p_thresh)
        group_vs_shuffle = bh_enrichment_df.loc[ind] < curr_thresh
        sig_pairs = group_vs_shuffle[
            group_vs_shuffle]  # .stack().reset_index().rename({"level_1":atac_col, 0:"p_val"}, axis=1)
        # sig_pairs = sig_pairs.sort_values("p_val")
        #     sig_pairs
        clone_min_sig[ind] = sig_pairs

        sns.histplot(data=data, stat='probability', color='blue',
                     ax=axs[curr_row, curr_col], binwidth=0.05, binrange=(0,1))

        sns.histplot(bh_enrichment_df.loc[ind], stat='probability', binwidth=0.05, binrange=(0,1),
                     color='green', alpha=0.05,
                     ax=axs[curr_row, curr_col])
        # sns.histplot(data=real_df.loc[data["variable"][0]], x='value',stat='probability', color='green', alpha=0.05)
        plt.xlabel("P-value")
        # print(real_df.loc[data["variable"][0]])
        # plt.title(ind)
        #print('curr_thresh', curr_thresh)
        axs[curr_row, curr_col].axvline(curr_thresh, color='r')

    f.tight_layout()
    if out_f is not None:
        plt.savefig(out_f,dpi=300, bbox_inches = "tight")

    return clone_min_sig, shuffle_pval_df, curr_thresh


def get_out(shuffle, clones, bh_enrichment_df, p_thresh, clone_map, atac_col, outdir, figs_close=True):
    global_min, clone_min, global_all, clone_all = get_shuffle_results(shuffle, clones, clone_map=clone_map)
    #ic('global all')
    out_all = plot_glob_all(bh_enrichment_df, global_all, p_thresh, atac_col=atac_col, out_f=join(outdir, "shuffle_all.pdf"))
    #ic('global min')
    out_min = plot_glob_min(bh_enrichment_df, global_min, p_thresh, atac_col=atac_col, out_f=join(outdir, "shuffle_min.pdf"))
    #ic('clone all')
    out_cloneall = plot_clone_all(bh_enrichment_df, clone_all, p_thresh, out_f=join(outdir, "shuffle_cloneMin.pdf"))
    #ic('clone min')
    out_clonemin = plot_clone_min(bh_enrichment_df, clone_min, p_thresh=p_thresh, out_f=join(outdir, "shuffle_cloneAll.pdf"))
    if figs_close:
        plt.close('all')
    # save results df with p-vals for each
    #results_df = pd.concat((out_all[1].melt(), out_min[1].melt(), out_cloneall[1].melt(), out_clonemin[1].melt()))
    results_df = pd.concat((out_all[1].reset_index().melt(id_vars=["index"]).assign(method="global_all"),out_min[1].reset_index().melt(id_vars=["index"]).assign(method="global_min"),
                            out_cloneall[1].reset_index().melt(id_vars=["index"]).assign(method="clone_all"),
                            out_clonemin[1].reset_index().melt(id_vars=["index"]).assign(method="clone_min")),
                           axis=0)

    out_d = {"sig_all":out_all, "sig_min":out_min,
             "sig_cloneAll": out_cloneall, "sig_cloneMin": out_clonemin}
    pickle.dump(out_d, open(join(outdir, "shuffle_results.p"), "wb"))
    results_df.to_csv(join(outdir, "shuffle_results_pvals.tsv"),sep="\t")
    return results_df, out_d



def plot_sig_results_per_method(results_df, p_thresh, outdir):
    for ind, val in results_df[results_df["value"] < p_thresh].groupby("method"):
        if len(val)==0 or val.shape[1]==0:
            print('val empty')
            print(ind)
            continue

        curr_df = val.astype(object).pivot(index="index", columns="variable",values="value").fillna(1)
        if (curr_df.shape[0] > 1) and (curr_df.shape[1] > 1):
            g = sns.clustermap(curr_df, vmax=float(p_thresh)+0.05)
        else:
            f = plt.figure()
            sns.heatmap(curr_df, vmax=float(p_thresh)+0.05)

        title = f"{ind}\np-val;cutoff is {p_thresh}; vmax {float(p_thresh)+0.05}"
        plt.suptitle(title)
        #g.fig.suptitle(title)
        #g.ax_cbar.set(title=title)  # g.ax_cbar.set(title="-log10 p-value")
        plt.savefig(join(outdir, f"{ind}_shuffle_sig"),dpi=300, bbox_inches = "tight")
    return


def hypergeo_plots(groups, clones, atac_cl, sizes, p_thresh, atac_col,
                   clone_col, outdir):
    print('hypergeo plots')
    print("Running hypergeo and saving sig results")
    print("plotting counts")
    groups["log2_count"] = np.log2(groups["count"] + 1)

    grp_counts = groups.pivot(index=atac_col, columns=clone_col, values="log2_count").fillna(0)
    print('grp_counts', grp_counts)
    if (len(groups[clone_col].unique()) > 1 and (grp_counts.shape[0] > 1)) and (grp_counts.shape[1] > 1):
        g = sns.clustermap(grp_counts)
    elif grp_counts.shape[0] > 0 and grp_counts.shape[1] > 0:
            f = plt.figure()
            sns.heatmap(grp_counts)
    else:
        plt.figure()
    plt.gca().set_title("log2 ncells")
    plt.savefig(join(outdir, "ncells.png"),dpi=300, bbox_inches = "tight")
    groups.rename({atac_col: "clusterID", clone_col: "cloneID"}, axis=1).to_csv(join(outdir, "ncells.csv"))

    output_df, bh_enrichment_df = pipeline_groups_hypergeo(groups,
                                                           clones,
                                                           atac_cl,
                                                           sizes,
                                                           p_thresh,
                                                           atac_col,
                                                           clone_col,
                                                           to_filt=False,
                                                           to_correct=True)
    bh_enrichment_df.to_csv(join(outdir, "hypergeo_padjusted.csv"))
    output_df.to_csv(join(outdir, "hypergeo_padjusted_sigOnly.csv"))

    if output_df.shape[0] == 0 and output_df.shape[1] == 0:
        f, ax = plt.subplots()
        plt.title("Empty data (none are significant)")
        f.savefig(join(outdir, "hypergeo_padjusted.png"), dpi=300,
                    bbox_inches="tight")
        return
    elif output_df.shape[0] <= 1 or (bh_enrichment_df.loc[output_df.index].shape[0] <= 1):
        plt.figure()
        sns.heatmap(-np.log10(bh_enrichment_df.fillna(1)))
        plt.title("-log10 BH-adjusted p-values for clonal shifts across clusters.")
        plt.savefig(join(outdir, "hypergeo_padjusted.png"), dpi=300,
                    bbox_inches="tight")
    else:
        g = sns.clustermap(
            -np.log10(bh_enrichment_df.loc[output_df.index].fillna(1)),
            row_cluster=False)
        g.ax_heatmap.set(xlabel="Cluster ID")
        g.ax_cbar.set(title="-log10 p-value")
        g.fig.suptitle("-log10 BH-adjusted p-values for clonal shifts across clusters.")
        plt.savefig(join(outdir, "hypergeo_padjusted.png"),dpi=300, bbox_inches = "tight")
    return True


def run_data_and_shuffle(groups, outdir, atac_col, clone_col, p_thresh, clones, atac_cl, to_p_correct=False, n_shuffle=1000, figs_close=False, n_cpus=4):
    # Run enrichment
    hyper_df = create_enrichment(groups, atac_col, clone_col, p_thresh, clones=clones, atac_cl=atac_cl, to_correct=False)

    # Run shuffled enrichment-same thing but shuffled labels and ran n_shuffle times
    shuffle = shuffle_hypergeo(groups, atac_col, clone_col, p_thresh, clones, atac_cl, to_p_correct=to_p_correct, n_shuffle=n_shuffle,
                                  to_parallel=True, n_cpus=n_cpus)

    # Get the output of the shuffle results and the associated p-values
    clone_map = {x: ind for ind, x in enumerate(clones)}
    results_df, out_d = get_out(shuffle, clones, hyper_df, p_thresh, clone_map,
                                atac_col, outdir, figs_close=figs_close)

    out_df = results_df.copy()
    out_df = out_df[out_df["value"] < p_thresh]
    print('out_df', out_df)
    # out_df["BH_p_adj"] = out_df.apply(
    #     lambda x: hyper_df.loc[(x["clone"], x["lineage"]), "BH_p_adj"],
    #     axis=1)
    out_df = out_df.sort_values(["method", "value"]) #, "BH_p_adj"])
    out_df.to_csv(join(outdir, "sig_results.tsv"))
    plot_sig_results_per_method(results_df, p_thresh, outdir)

    return out_df, hyper_df, results_df, out_d


def check_sig(x):
    #print(x.head())
    print('name', x.name)
    name, clust, cond = x.name
    x = x.set_index("method")["is_sig"]

    assert(x.index.duplicated().sum()==0)

    sig = 0
    if 'hypergeo' not in x.index:
        print('hypergeo not sig')
        return sig
    if x["hypergeo"]==True:
        if ("global_all" in x) and x["global_all"]==True:
            if ("clone_min" in x) and (x["clone_min"]==True):
                if ("global_min" in x) and (x["global_min"]==True):
                    sig = 4
                else:
                    sig = 3
            else:
                sig = 2
        else:
            sig = 1
    return sig


def hypergeo_score(df, p_thresh):
    p_df = df.copy()
    p_df["is_sig"] = p_df["pval"]<p_thresh
    print(p_df.shape)
    print(p_df.duplicated(["index","variable", "condition", "method"]).any())

    p_df_group = p_df.groupby(["index", "variable", "condition"]).apply(
        check_sig)
    print(p_df_group.shape)
    print(p_df_group)
    print(p_df_group.index.name)
    p_df_out = p_df_group
    if "index" not in p_df_out:
        p_df_out = p_df_out.reset_index()

    # p_df_out = p_df_group.reset_index().rename({0:"significant_score", "index":"clone", "variable":"cluster"}, axis=1)
    print('p_df_out')
    print(p_df_out.head())
    if "significant_score" not in p_df_out:
        p_df_out = p_df_out.rename({0: "significant_score"},
                                               axis=1)
    print('p_df_out after')
    print(p_df_out.head())
    p_df_out["cluster_condition"] = p_df_out.apply(
        lambda x: f'{x["variable"]}_{x["condition"]}', axis=1)

    p_df_out = p_df_out.pivot(index='index',
                              columns='cluster_condition',
                              values="significant_score").fillna(0)

    return p_df_group, p_df_out


