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


def run_enrichment(df):
    """Runs hypergeometric for flt3 expansion.

    df: pd.DataFrame where index is Control and Flt3, columns are cluster labels, and elements are number of cells.
    """
    enrichment_df = pd.DataFrame(index=["p", "Fisher p"],
                                 columns=df.columns, dtype=np.float128)

    # M: Total number of cells
    # n: Number of cells in the clone
    # N: Number of flt3 cells
    # x: Number of cells in specific clone with flt3
    M = df.sum().sum()
    N = df.loc["# Flt3 Cells in Cluster"].sum()  # df.loc["Flt3"].sum()
    # rv = hypergeom(M, n, N)
    for col in df.columns:
        n = df[col].sum()
        # x = df.loc["Flt3", col]
        x = df.loc["# Flt3 Cells in Cluster", col]
        prb = 1 - hypergeom.cdf(x, M, n, N)
        enrichment_df.loc["p", col] = prb

        oddsratio, fish_p = fisher_exact(
            [[x, N - x], [n - x, M - N - (n - x)]])
        enrichment_df.loc["Fisher p", col] = fish_p

    return enrichment_df

@click.command()
@click.argument("donor_indir",  type=click.Path(exists=True))
@click.argument("clones_indir",  type=click.Path(exists=True))
@click.argument("outdir",  type=click.Path(exists=True))
@click.argument("n_donors",type=click.INT)
@click.argument("n_clone_list",type=click.INT, nargs=-1)
def main(donor_indir, clones_indir, outdir, n_donors, n_clone_list):
    for n in range(n_donors):
        print("Donor N")
        curr_ad_f = join(donor_indir, f"donor{n}.AD.mtx")
        curr_dp_f = join(donor_indir, f"donor{n}.DP.mtx")
        print(curr_ad_f)
        print(curr_dp_f)
        print(n_clone_list)
        curr_ad = mmread(curr_ad_f).tocsc()
        curr_dp = mmread(curr_dp_f).tocsc()
        curr_labels = pd.read_csv(join(donor_indir, f"donor{n}.labels.txt"),
                                  index_col=0)['sample ID']
        print(set(curr_labels))
        #donor4.clones10.modelCA.p
        in_f= join(clones_indir, f"donor{n}")

        for k in n_clone_list:
            curr_modelCA = pickle.load(
                open(f"{in_f}.clones{k}.modelCA.p", "rb"))
            cell_clusters, sample_labels = extract_clusters_matrix(curr_modelCA,
                                                            curr_ad,
                                                            curr_dp,
                                                            curr_labels)
            # Create counts df
            clust_counts = pd.DataFrame(index=["# Control Cells in Cluster",
                                               "# Flt3 Cells in Cluster"],
                                        columns=sample_labels.keys())

            for curr_k in clust_counts.columns:
                clust_counts.at["# Control Cells in Cluster", curr_k] = (
                            sample_labels[curr_k][
                                "sample ID"] == "Control").sum()
                clust_counts.at["# Flt3 Cells in Cluster", curr_k] = (
                            sample_labels[curr_k][
                                "sample ID"] == "Flt3").sum()
            print(clust_counts )
            clust_counts = clust_counts.astype(np.double)
            fold_df = pd.DataFrame(np.log2(
                (clust_counts.loc["# Flt3 Cells in Cluster"] + 1) / (
                            clust_counts.loc[
                                "# Control Cells in Cluster"] + 1))).transpose()
            fold_df = fold_df.rename({0: "Flt3 fold enrichment"}, axis=0)

            # Get enrichment
            enrich_df = run_enrichment(clust_counts)
            enrich_stats = pd.concat(
                (clust_counts, fold_df, enrich_df)).transpose()

            f, ax = plt.subplots(nrows=3, ncols=1, figsize=(15, 15),
                                 dpi=300)
            labels = np.array(
                [f"{x:.2E}" for x in enrich_df.astype(float).values[0]])
            sns.heatmap(enrich_df.astype(float), annot=True, fmt=".2E",
                        cbar=False,  # labels.reshape([1, len(labels)]),
                        cmap='RdYlGn', ax=ax[0])

            ax[0].set_title("Flt3 enrichment p-value")

            sns.heatmap(fold_df, ax=ax[1], annot=True, fmt=".2f",
                        cbar=False, cmap="RdBu")
            ax[1].set_aspect('equal', adjustable='box')
            ax[1].set_title("Log2 Fold change of flt3/control")

            plt.xlabel("log2((Flt3+1)/(WT+1))")
            sns.heatmap((clust_counts.astype(int)), vmin=0, ax=ax[2],
                        annot=True, fmt=".1E", cbar=False, cmap="Blues")
            ax[2].set_title("Number of cells in each cluster")

            plt.suptitle(f"Donor {n}")
            plt.savefig(join(outdir, f"donor{n}.clones{k}_labelEnrich.png"))
            plt.close()

            # Put in scatterplot.
            enrich_stats = pd.concat(
                (clust_counts, fold_df, enrich_df)).transpose()
            enrich_stats['p'] = enrich_stats['p'].astype(np.float32)
            enrich_stats['-log10p'] = -np.log10(enrich_stats['p'])
            enrich_stats.loc[
                enrich_stats["-log10p"] == np.infty, "-log10p"] = 0
            f = plt.figure(figsize=(10, 10))
            sns.scatterplot(data=enrich_stats,
                            x="# Control Cells in Cluster",
                            y="# Flt3 Cells in Cluster", size="-log10p",
                            s=100, sizes=(20, 200),
                            hue="Flt3 fold enrichment", palette="RdBu")
            plt.savefig(join(outdir, f"donor{n}.clones{k}_scatterEnrich.png"))
            plt.close()
            enrich_stats['Fisher p'] = enrich_stats['Fisher p'].astype(
                np.float32)
            enrich_stats['Fisher -log10p'] = -np.log10(
                enrich_stats['Fisher p'])
            enrich_stats.loc[enrich_stats[
                                 "Fisher -log10p"] == np.infty, "Fisher -log10p"] = 0
            f, ax = plt.subplots(figsize=(10, 10))
            sns.scatterplot(data=enrich_stats,
                            x="# Control Cells in Cluster",
                            y="# Flt3 Cells in Cluster",
                            size="Fisher -log10p", s=100, sizes=(20, 200),
                            hue="Flt3 fold enrichment", palette="RdBu")
            plt.axis('square')
            ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='black')
            # ax.plot([0, max(max(enrich_stats[["Control", "Flt3"]]))])

            # helper_save(out_f+f"clones{k}_scatterFisherEnrich.png", to_pdf=False)
            plt.savefig(join(outdir, f"donor{n}.clones{k}_scatterFisherEnrich.png"),dpi=300)
            plt.close()

            enrich_stats.to_csv(join(outdir, "enrichment.csv"))
            # dummy output variable
            with open(join(outdir,".status"), 'w') as f:
                f.write('Completed')
            return


if __name__=="__main__":
    main()