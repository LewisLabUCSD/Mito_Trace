from os.path import join, exists, dirname
from glob import glob
import pickle
import mplh.cluster_help as ch
import os
#import vireoSNP
import numpy as np
from scipy import sparse
from scipy.io import mmread
import matplotlib.pyplot as plt
from scipy.stats import hypergeom, fisher_exact
#print(vireoSNP.__version__)
import click
import pandas as pd
import seaborn as sns
#from vireoSNP import Vireo
import matplotlib as mpl
mpl.use('Agg')

from icecream import ic
#np.set_printoptions(formatter={'float': lambda x: format(x, '.5f')})

#n_clone_list = [3, 5, 10, 20, 40]  # [2,3,4,5,6,7]

def extract_clusters_matrix(modelCA, ad, dp, sample_colors, prob_thresh=0.9,
                     doublet_thresh=0.9):
    """ Creates a dictionary where the keys are the cluster IDs (0-based index)
        and the values are the cell indices (0-based index)"""
    cell_clusters = {}
    sample_labels = {}
    doublet_prob = \
    modelCA.predict_doublet(ad, dp, update_GT=False, update_ID=False)[
        0].sum(axis=1)
    low_conf_cells = np.flatnonzero(doublet_prob > doublet_thresh)
    for n in range(modelCA.ID_prob.shape[1]):
        cell_clusters[n] = np.flatnonzero(
            (modelCA.ID_prob[:, n] > prob_thresh))
        print('before doublet removal')
        print(len(cell_clusters[n]))
        cell_clusters[n] = cell_clusters[n][
            ~(np.isin(cell_clusters[n], low_conf_cells))]
        print('after doublet removal')
        print(len(cell_clusters[n]))
        curr_sample_colors = sample_colors.iloc[cell_clusters[n]].copy()
        curr_sample_colors = curr_sample_colors.reset_index()
        sample_labels[n] = curr_sample_colors
        print(f"Cluster {n}: {len(cell_clusters[n])} cells ")

    return cell_clusters, sample_labels


def run_enrichment_stats(df, flt_var="Flt3l"):
    """Runs hypergeometric for flt3 expansion.

    df: pd.DataFrame where index is Control and Flt3, columns are cluster labels, and elements are number of cells.
    """
    enrichment_df = pd.DataFrame(index=["hypergeom p", "Fisher p"],
                                 columns=df.columns, dtype=np.float128)
    ic('df')
    ic(df.head())
    # M: Total number of cells
    # n: Number of cells in the clone
    # N: Number of flt3 cells
    # x: Number of cells in specific clone with flt3
    M = df.sum().sum()
    N = df.loc[f"# {flt_var} Cells in Cluster"].sum()  # df.loc["Flt3"].sum()
    # rv = hypergeom(M, n, N)
    for col in df.columns:
        n = df[col].sum()
        # x = df.loc["Flt3", col]
        x = df.loc[f"# {flt_var} Cells in Cluster", col]
        prb = 1 - hypergeom.cdf(x, M, n, N)
        enrichment_df.loc["hypergeom p", col] = prb

        oddsratio, fish_p = fisher_exact(
            [[x, N - x], [n - x, M - N - (n - x)]])
        enrichment_df.loc["Fisher p", col] = fish_p
    enrichment_df=enrichment_df.transpose()
    enrichment_df['Fisher p'] = enrichment_df['Fisher p'].astype(
        np.float32)
    enrichment_df['Fisher -log10p'] = -np.log10(enrichment_df['Fisher p'])
    enrichment_df.loc[enrichment_df["Fisher -log10p"] == np.infty, "Fisher -log10p"] = 0
    return enrichment_df.transpose()


def wrap_lineage_enrichment(clones_indir, outdir, n_clone_list, samples,
                       plot_ind, name="", names=None):
    for nclones in n_clone_list:
        lineage_enrichment(clones_indir, outdir, nclones, samples,
                       plot_ind, name="", names=None)



def lineage_enrichment(clones_indir, outdir, nclones, samples,
                       plot_ind, name="", names=None, pseudocount=1):

    if names is None:
        names = samples
    #all_enrich = {}
    all_enrich_norm ={}

    ###############################################################
    ## Added 06/08
    ###############################################################
    cells_meta = pd.read_csv(join(clones_indir, "cells_meta.tsv"),
                             sep='\t')
    cells_meta["lineage"] = cells_meta["lineage"].astype('Int64')
    cells_meta["donor"] = cells_meta["donor"].astype('Int64')
    for d, curr_donor in cells_meta.groupby("donor"):
        # Create counts df
        print('d', d)
        print(curr_donor.groupby(["condition", "lineage"]).size())

        clust_counts = curr_donor.groupby(["condition", "lineage"]).size().reset_index().pivot(index='condition',columns='lineage', values=0).fillna(0)
        clust_counts = clust_counts + pseudocount

        if len(clust_counts) == 0:
            print("No lineages detected in donor. Continuing")
            continue
        clust_counts = clust_counts.rename({x:y for (x,y) in zip(samples, names)}, axis=0)

        #clust_counts.columns = clust_counts.columns.astype('Int64')
        clust_counts.index = [f"# {x} Cells in Cluster" for x in clust_counts.index]
        ###############################################################
        clust_counts = clust_counts.astype('Int64')
        print('clust_counts')
        print(clust_counts)

        # Get enrichment
        if "Input" in names:
            enrich_df = run_enrichment_stats(clust_counts.drop("# Input Cells in Cluster", axis=0), flt_var=names[1])
        else:
            enrich_df = run_enrichment_stats(clust_counts, flt_var=names[1])

        ### Convert cluster numbers into probablilites.
        clust_counts_norm = (clust_counts).copy().astype(np.double)
        clust_counts_norm  = clust_counts_norm.div(clust_counts_norm.sum(axis=1), axis='rows')
        # print('clust_counts_norm')
        # print(clust_counts_norm)
        fold_df_norm = pd.DataFrame(
                    (clust_counts_norm.loc[f"# {names[1]} Cells in Cluster"]) / (
                        clust_counts_norm.loc[f"# {names[0]} Cells in Cluster"])).transpose()
        fold_df_norm = fold_df_norm.rename({0: f"{names[1]} fold enrichment norm"}, axis=0)

        print('fold_df_norm')
        print(fold_df_norm )
        enrich_stats = create_enrich(clust_counts, fold_df_norm,  enrich_df)
        enrich_stats.to_csv(join(outdir, f"{name}enrichmentNorm_clones{nclones}_donor{d}.csv"))
        if plot_ind:
            plot_volcano(enrich_stats, f_save=join(outdir, f"{name}volcano_donor{d}.clones{nclones}_Fisher_foldNorm.png"), v=1, x=f"{names[1]} fold enrichment norm")
            plot_volcano(enrich_stats, f_save=join(outdir,f"{name}volcano_donor{d}.clones{nclones}_hyper_foldNorm.png"),y="-log10p", v=1, x=f"{names[1]} fold enrichment norm")
        all_enrich_norm[(d, nclones)] = enrich_stats

    all_enrich_df = pd.concat(all_enrich_norm).reset_index().rename(
        {"level_0": "Donor", "level_1": "clones"}, axis=1)

    #print("Does donor have any cells")
    #print("Any null"(all_enrich_df["Donor"].isnull()).any())
    all_enrich_df["Donor"] = all_enrich_df["Donor"].astype(object)

    for nclones_ind, val in all_enrich_df.groupby("clones"):
        # Main figure is Fisher norm results
        plot_volcano(val, hue="Donor", x=f"{names[1]} fold enrichment norm",
                     size=f"# {names[1]} Cells in Cluster",
                     f_save=join(outdir, f"{name}volcano_Fisher_foldNorm.png"), v=1)
        # # Also plot without the norm
        # plot_volcano(val, hue="Donor", size=f"# {names[1]} Cells in Cluster",
        #              f_save=join(outdir, f"{name}volcano_Fisher_fold.png"), v=1, x=f"{names[1]} fold enrichment")
        # Plot the hypergeometric p-value
        plot_volcano(val, hue="Donor", x=f"{names[1]} fold enrichment norm",
                     size=f"# {names[1]} Cells in Cluster",
                     f_save=join(outdir, f"{name}volcano_HyperG_foldNorm.png"), v=1, y = "-log10p")
        # Have the size be based on number of cells in baseline sample
        plot_volcano(val, hue="Donor", x=f"{names[1]} fold enrichment norm",
                     size=f"# {names[0]} Cells in Cluster",
                     f_save=join(outdir, f"{name}volcano_Fisher_foldNorm_{names[0]}Size.png"), v=1)


    # dummy output variable
    with open(join(outdir,".status"), 'w') as f:
        f.write('Completed')
    return


def plot_volcano(enrich_stats, x="Flt3l fold enrichment",
                 y="Fisher -log10p", hue=None, f_save=None, v=0,
                 size=None, to_close=True, to_log=False, ylim=None,
                 xlim=None):
    enrich_stats = enrich_stats.astype(float)

    f, ax = plt.subplots(figsize=(10, 10), sharey=True, sharex=True)

    if to_log:
        enrich_stats[f"log2 {x}"] = np.log2(enrich_stats[x])
        x = f"log2 {x}"
        v=0

    if hue is None:
        n_clr=1
    else:
        n_clr=len(set(enrich_stats[hue]))
    sns.scatterplot(data=enrich_stats, x=x,
                    y=y, s=100, sizes=(20, 200),
                    palette=sns.color_palette("Set1", n_clr),
                    ax=ax, hue=hue, size=size)
    plt.axvline(x=v)
    #ax.plot([0.5], [0.5], transform=ax.transAxes, color='black')
    #plt.axis('square')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
               borderaxespad=0)
    if f_save is not None:
        plt.savefig(f_save, bbox_inches='tight')
    if to_close:
        plt.close()
    return


def wrap_plot_volcano(enrich_stats_all, x="Flt3l fold enrichment norm",
                 y="Fisher -log10p", hue=None, f_save=None, v=0,
                 size=None, to_close=False, to_log=False):

    f, ax = plt.subplots(figsize=(15, 15), ncols=len(enrich_stats_all),
                         sharey=True, sharex=True)

    count = 0
    for k in enrich_stats_all:
        curr_ax = ax[count]

        enrich_stats = enrich_stats_all[k].copy()
        enrich_stats = enrich_stats.astype(float)

        if to_log:
            new_x = f"log2 {x}"
            enrich_stats[new_x] = np.log2(enrich_stats[x])
            v = 0
        else:
            new_x = x

        if hue is None:
            n_clr = 1
        else:
            n_clr = len(set(enrich_stats[hue]))
        sns.scatterplot(data=enrich_stats, x=new_x, y=y, s=100, sizes=(20, 200),
                        palette=sns.color_palette("Set1", n_clr),
                        ax=curr_ax, hue=hue, size=size)

        curr_ax.axvline(x=v)
        # ax.plot([0.5], [0.5], transform=ax.transAxes, color='black')
        # plt.axis('square')
        curr_ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
                  borderaxespad=0)
        curr_ax.set_title(k)
        count+=1
    if f_save is not None:
        plt.tight_layout()
        plt.savefig(f_save, bbox_inches='tight', dpi=300)
    if to_close:
        plt.close()
    return


def create_enrich(clust, fold, enrich, fold_norm=None):
    # print('clust')
    # print(clust.head())
    # print('fold')
    # print(fold.head())
    # print('enrich')
    # print(enrich.head())
    if fold_norm is not None:
        enrich_stats = pd.concat((clust.transpose(), fold.transpose(),
                                  enrich.transpose(), fold_norm.transpose()),
                                  axis=1)
    else:
        enrich_stats = pd.concat((clust.transpose(), fold.transpose(),
                                  enrich.transpose()),
                                  axis=1)
        # enrich_stats = pd.concat(
        #     (clust, fold, enrich)).transpose()
    enrich_stats['hypergeom p'] = enrich_stats['hypergeom p'].astype(float)
    enrich_stats['-log10p'] = -np.log10(enrich_stats['hypergeom p'])
    enrich_stats.loc[enrich_stats["-log10p"] == np.infty, "-log10p"] = 0
    return enrich_stats


def plot_label_enrich(enrich_df, fold_df, clust_counts, outdir, n, k):
    f, ax = plt.subplots(nrows=3, ncols=1, figsize=(15, 15), dpi=300)
    labels = np.array(
        [f"{x:.2E}" for x in enrich_df.astype(float).values[0]])
    sns.heatmap(enrich_df.astype(float), annot=True, fmt=".2E",
                cbar=False,  # labels.reshape([1, len(labels)]),
                cmap='RdYlGn', ax=ax[0])
    ax[0].set_title("Flt3 enrichment p-value")
    #print('fold_df', fold_df)
    sns.heatmap(fold_df, ax=ax[1], annot=True, fmt=".2f", cbar=False,
                cmap="RdBu")
    ax[1].set_aspect('equal', adjustable='box')
    ax[1].set_title("Fold change of flt3/control")

    plt.xlabel("((Flt3+1)/(WT+1))")
    sns.heatmap((clust_counts.astype(int)), vmin=0, ax=ax[2],
                annot=True, fmt=".1E", cbar=False, cmap="Blues")
    ax[2].set_title("Number of cells in each cluster")
    plt.suptitle(f"Donor {n}")
    plt.savefig(join(outdir, f"donor{n}.clones{k}_labelEnrich.png"))
    plt.close()

    return


def plot_cluster_nums(enrich_stats, f_save, samples):
    f, ax = plt.subplots(figsize=(10, 10))
    sns.scatterplot(data=enrich_stats,
                    x=f"# {samples[0]} Cells in Cluster",
                    y=f"# {samples[1]} Cells in Cluster",
                    size="Fisher -log10p", s=100, sizes=(20, 200),
                    hue=f"{samples[1]} fold enrichment", palette="RdBu")
    plt.axis('square')
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='black')
    plt.savefig(f_save, dpi=300)
    # plt.savefig(join(outdir, f"donor{n}.clones{k}_scatterFisherEnrich.png"),dpi=300)
    plt.close()
    return


def plot_clone_scatter(enrich_stats, outdir, samples, n, k, size='-log10p'):
    f = plt.figure(figsize=(10, 10))
    sns.scatterplot(data=enrich_stats,
                    x=f"# {samples[0]} Cells in Cluster",
                    y=f"# {samples[1]} Cells in Cluster",
                    size=size, s=100, sizes=(20, 200),
                    hue=f"{samples[1]} fold enrichment", palette="RdBu")
    #plt.show()
    plt.savefig(join(outdir, f"donor{n}.clones{k}_scatterEnrich.png"))
    plt.close()


@click.command()
@click.argument("clones_indir", type=click.Path(exists=True))
@click.argument("outdir",  type=click.Path(exists=True))
@click.argument("nclones",type=click.INT)
@click.argument("samples", type=click.STRING)
@click.option("--plot_ind", default=False)
@click.option("--tests", default="")
def main(clones_indir, outdir, nclones, samples, plot_ind, tests):
    if tests != "":
        comps = tests.split(";")
        for c in comps:
            name, a_name, a, b_name, b = c.split(',')
            samples = [a, b]
            lineage_enrichment(clones_indir, outdir, nclones,
                               samples, plot_ind, name=f"{name}_",
                               names=[a_name, b_name])
    else:
        samples = samples.split(',')
        lineage_enrichment(clones_indir, outdir, nclones, samples, plot_ind)
    return


if __name__=="__main__":
    main()
