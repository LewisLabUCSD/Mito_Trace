from os.path import join, exists, dirname
from glob import glob
import pickle
import mplh.cluster_help as ch
import os
import vireoSNP
import numpy as np
from scipy import sparse
from scipy.io import mmread
import matplotlib.pyplot as plt
from scipy.stats import hypergeom, fisher_exact
print(vireoSNP.__version__)
import click
import pandas as pd
import seaborn as sns
from vireoSNP import Vireo
import matplotlib as mpl
mpl.use('Agg')
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


def run_enrichment(df, flt_var="Flt3"):
    """Runs hypergeometric for flt3 expansion.

    df: pd.DataFrame where index is Control and Flt3, columns are cluster labels, and elements are number of cells.
    """
    enrichment_df = pd.DataFrame(index=["hypergeom p", "Fisher p"],
                                 columns=df.columns, dtype=np.float128)

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



def lineage_enrichment(clones_indir, outdir, n_clone_list, samples, plot_ind):
    #n_clone_list=[20,100]
    samples = samples.split(",")
    all_enrich = {}
    all_enrich_norm ={}

    ###############################################################
    ## Added 06/08
    ###############################################################
    for k in n_clone_list:
        cells_meta = pd.read_csv(join(clones_indir, f"lineage{k}", "cells_meta.tsv"),
                                 sep='\t')
        for n, curr_donor in cells_meta.groupby("donor"):
            # Create counts df
            clust_counts = curr_donor.groupby(["condition", "lineage"]).size().reset_index().pivot(index='condition',columns='lineage', values=0).fillna(0)
            clust_counts.index = [f"# {x} Cells in Cluster" for x in clust_counts.index]

    ###############################################################
    ###############################################################
    ## Replaces this part
    # for n in range(n_donors):
    #     print("Donor N")
    #     curr_ad_f = join(donor_indir, f"donor{n}.AD.mtx")
    #     curr_dp_f = join(donor_indir, f"donor{n}.DP.mtx")
    #
    #     curr_ad = mmread(curr_ad_f).tocsc()
    #     curr_dp = mmread(curr_dp_f).tocsc()
    #     curr_labels = pd.read_csv(join(donor_indir, f"donor{n}.labels.txt"),
    #                               index_col=0)['sample ID']
    #     #print(set(curr_labels))
    #     #donor4.clones10.modelCA.p
    #
    #     for k in n_clone_list:
    #         in_f = join(clones_indir, f"donor{n}")
    #         curr_modelCA = pickle.load(
    #             open(f"{in_f}.clones{k}.modelCA.p", "rb"))
    #
    #         cell_clusters, sample_labels = extract_clusters_matrix(curr_modelCA,
    #                                                         curr_ad,
    #                                                         curr_dp,
    #                                                         curr_labels)
    #
    #         # Create counts df
    #         index = [f"# {x} Cells in Cluster" for x in samples]
    #         clust_counts = pd.DataFrame(index=index,
    #                                     columns=sample_labels.keys())
    #         for curr_k in clust_counts.columns:
    #             for s in samples:
    #                 clust_counts.at[
    #                     f"# {s} Cells in Cluster", curr_k] = (
    #                         sample_labels[curr_k]["sample ID"] == s).sum()
            ###############################################################
            clust_counts = clust_counts.astype('Int64')
            # Get enrichment
            #print('clust_conts')
            #print(clust_counts.index)
            if "Input" in samples:
                enrich_df = run_enrichment(clust_counts.drop("# Input Cells in Cluster", axis=0), flt_var=samples[1])
            else:
                enrich_df = run_enrichment(clust_counts, flt_var=samples[1])
            print('enrich_df')
            print(enrich_df)
            fold_df = pd.DataFrame(
                    (clust_counts.loc[f"# {samples[1]} Cells in Cluster"] + 1) / (
                            clust_counts.loc[f"# {samples[0]} Cells in Cluster"] + 1)).transpose()
            fold_df = fold_df.rename({0: "Flt3l fold enrichment"}, axis=0)
            plot_label_enrich(enrich_df, fold_df, clust_counts, outdir, n, k)
            enrich_stats = create_enrich(clust_counts, fold_df, enrich_df)
            enrich_stats.to_csv(join(outdir, "enrichment.csv"))
            all_enrich[(n, k)] = enrich_stats
            if plot_ind:
                plot_volcano(enrich_stats, f_save=join(outdir,f"volcano_donor{n}.clones{k}_Fisher_fold.png"))
                plot_volcano(enrich_stats, f_save=join(outdir,f"volcano_donor{n}.clones{k}_hyper_fold.png"), y="Fisher -log10p")

            ### Convert cluster numbers into probablilites.
            #print('clust_counts', clust_counts)
            clust_counts_norm = (clust_counts+1).copy().astype(np.double)
            clust_counts_norm  = clust_counts_norm.div(clust_counts_norm.sum(axis=1), axis='rows')
            print('clust_counts_norm')
            print(clust_counts_norm)
            fold_df_norm = pd.DataFrame(
                        (clust_counts_norm.loc[f"# {samples[1]} Cells in Cluster"] + 1) / (
                            clust_counts_norm.loc[f"# {samples[0]} Cells in Cluster"] + 1)).transpose()
            fold_df_norm = fold_df_norm.rename({0: "Flt3l fold enrichment norm"}, axis=0)

            #fold_df_norm = (clust_counts_norm).transpose()
            #fold_df_norm = fold_df_norm.rename({0: "Flt3l fold enrichment norm"}, axis=1).transpose()
            #fold_df_norm = fold_df_norm.rename({0: "Flt3l fold enrichment norm"}, axis=1).transpose()
            print('fold_df_norm')
            print(fold_df_norm )
            enrich_stats = create_enrich(clust_counts, fold_df_norm,  enrich_df, fold_df)

            #enrich_stats = pd.concat((enrich_stats.transpose(), fold_df)).transpose()
            enrich_stats.to_csv(join(outdir, "enrichmentNorm.csv"))
            print('enrich_stats')
            print(enrich_stats)
            print(enrich_stats.columns)
            if plot_ind:
                plot_volcano(enrich_stats, f_save=join(outdir, f"volcano_donor{n}.clones{k}_Fisher_foldNorm.png"), v=1)
                plot_volcano(enrich_stats, f_save=join(outdir,f"volcano_donor{n}.clones{k}_hyper_foldNorm.png"),y="-log10p", v=1)
            all_enrich_norm[(n, k)] = enrich_stats

    #all_enrich = pd.concat(all_enrich).reset_index()
    all_enrich_df = pd.concat(all_enrich_norm).reset_index().rename(
        {"level_0": "Donor", "level_1": "clones"}, axis=1)
    print((all_enrich_df["Donor"].isnull()).any())
    all_enrich_df["Donor"] = all_enrich_df["Donor"].astype(object)
    #print('all_enrich_df')
    #print(all_enrich_df.head())
    for k, val in all_enrich_df.groupby("clones"):
        print('k', k)
        print('val', val)
        plot_volcano(val, hue="Donor", size=f"# {samples[1]} Cells in Cluster",
                     f_save=join(outdir, f"volcano.clones{k}_Fisher_fold.png"), v=1)
        plot_volcano(val, hue="Donor", size=f"# {samples[1]} Cells in Cluster",
                     f_save=join(outdir, f"volcano.clones{k}_HyperG_fold.png"), v=1, y = "-log10p")

        plot_volcano(val, hue="Donor", x="Flt3l fold enrichment norm",
                     size=f"# {samples[1]} Cells in Cluster",
                     f_save=join(outdir, f"volcano.clones{k}_Fisher_foldNorm.png"), v=1)
        plot_volcano(val, hue="Donor", x="Flt3l fold enrichment norm",
                     size=f"# {samples[1]} Cells in Cluster",
                     f_save=join(outdir, f"volcano.clones{k}_HyperG_foldNorm.png"), v=1, y = "-log10p")
        if "# Input Cells in Cluster" in val.columns:
            plot_volcano(val, hue="Donor", x="Flt3l fold enrichment norm",
                         size="# Input Cells in Cluster",
                         f_save=join(outdir, f"volcano.clones{k}_Fisher_foldNorm_inputSize.png"), v=1)
            plot_volcano(val, hue="Donor", x="Flt3l fold enrichment norm",
                         size="# Input Cells in Cluster",
                         f_save=join(outdir, f"volcano.clones{k}_HyperG_foldNorm_inputSize.png"), v=1, y = "-log10p")
        plot_volcano(val, hue="Donor", x="Flt3l fold enrichment norm",
                     size=f"# {samples[0]} Cells in Cluster",
                     f_save=join(outdir, f"volcano.clones{k}_Fisher_foldNorm_ControlSize.png"), v=1)
    # dummy output variable
    with open(join(outdir,".status"), 'w') as f:
        f.write('Completed')
    return


@click.command()
@click.argument("clones_indir", type=click.Path(exists=True))
@click.argument("outdir",  type=click.Path(exists=True))
@click.argument("n_clone_list",type=click.INT, nargs=-1)
@click.argument("samples", type=click.STRING)
@click.option("--plot_ind", default=True)
def main(clones_indir, outdir, n_clone_list, samples, plot_ind):
    lineage_enrichment(clones_indir, outdir, n_clone_list, samples, plot_ind)

    return


def plot_volcano(enrich_stats, x="Flt3l fold enrichment",
                 y="Fisher -log10p", hue=None, f_save=None, v=0,
                 size=None):
    f, ax = plt.subplots(figsize=(10, 10))
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
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
               borderaxespad=0)
    if f_save is not None:
        plt.savefig(f_save)#, bbox_inches='tight')
    plt.close()
    return


def create_enrich(clust, fold, enrich, fold_norm=None):
    if fold_norm is not None:
        enrich_stats = pd.concat(
            (clust, fold, enrich, fold_norm)).transpose()
    else:
        enrich_stats = pd.concat(
            (clust, fold, enrich)).transpose()
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

    # print('enrich_stats')
    # print(enrich_stats.head())

    return


def plot_cluster_nums(enrich_stats, f_save, samples):
    f, ax = plt.subplots(figsize=(10, 10))
    sns.scatterplot(data=enrich_stats,
                    x=f"# {samples[0]} Cells in Cluster",
                    y=f"# {samples[1]} Cells in Cluster",
                    size="Fisher -log10p", s=100, sizes=(20, 200),
                    hue="Flt3l fold enrichment", palette="RdBu")
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
                    hue="Flt3l fold enrichment", palette="RdBu")
    #plt.show()
    plt.savefig(join(outdir, f"donor{n}.clones{k}_scatterEnrich.png"))
    plt.close()




if __name__=="__main__":
    main()