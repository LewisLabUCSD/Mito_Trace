#import pymc3 as pm
import numpy as np
from numpy import random
import os
import pandas as pd
from tqdm import tqdm
#from src.config import ROOT_DIR
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import glob
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, roc_curve, average_precision_score, confusion_matrix
from scipy.spatial.distance import cdist
from pandarallel import pandarallel


from mplh.color_utils import get_colors
from mplh.fig_utils import legend_from_color
from mplh import cluster_help as ch
from src.simulations.utils.config import read_config_file, write_config_file

from dynamicTreeCut import cutreeHybrid
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage
from scipy.stats import entropy
from sklearn.model_selection import ParameterGrid


class Analysis:
    """Analysis of lineage tracing experiments

    Can look at both supervised (from simulation)- where the clones are
    known, along with unsupervised, which is the experimental side.
    """

    def __init__(self):
        return


    def calculate_average_clone_mt(self):
        return


    @staticmethod
    def cluster(cell_af):
        """Dynamic tree clustering of the rows of cell_af :param cell_af:
        :return:

        Args:
            cell_af:
        """
        distances = pdist(cell_af, "euclidean")
        link = linkage(distances, "average")
        clusters = cutreeHybrid(link, distances)['labels']
        return clusters

    @staticmethod
    def cluster_kmeans(X, min_c=2, max_c=10, n_clust=None):
        """
        Args:
            cell_af:
        """

        best_n_clusters = -1
        sil_score_max = -1  # this is the minimum possible score

        if n_clust is None:
            for n_clusters in range(min_c, max_c):
                model = KMeans(n_clusters=n_clusters, init='k-means++',
                               max_iter=100, n_init=1)
                labels = model.fit_predict(X)
                sil_score = silhouette_score(X, labels)
                print(
                    "The average silhouette score for %i clusters is %0.2f" % (
                    n_clusters, sil_score))
                if sil_score > sil_score_max:
                    sil_score_max = sil_score
                    best_n_clusters = n_clusters

            assert(best_n_clusters >= 0)

        else:
            best_n_clusters = n_clust

        model = KMeans(n_clusters=best_n_clusters, init='k-means++',
                       max_iter=100, n_init=1)
        labels = model.fit_predict(X)

        return labels


    @staticmethod
    def cluster_performance(y_true, y_pred):
        """

        :param true: vector of true clone labels
        :type true: iterable
        :param pred: vector of predicted clone labels
        :type pred: iterable
        :return: metrics dataframe for cluster performance
        :rtype: pd.DataFrame
        """
        metrics = dict()
        dropout = []
        rocs = []
        prec_scores = []
        prec_scores.append(average_precision_score(y_true, y_pred))
        return pd.DataFrame(metrics)


    @staticmethod
    def estimate_growth_rate(cells_meta, clone_col='clone',
                             pseudo_count=1):
        """ Estimates the growth rate for each clone, which in this case
        is based on the prediction

        :param cells_meta: Information on cells
            pd.DataFrame:
                Index:
                    RangeIndex
                Columns:
                    * Name: clone, dtype: object
                    * Name: pre_post, dtype: bool
        :return: clone_growth
        :rtype: pd.Series
            Name: RangeIndex, dtype: float
        """

        clone_sizes = cells_meta.groupby([clone_col, 'pre_post']).size().rename('Count').reset_index()
        #clone_sizes["Count"] += pseudo_count
        growth_estimate = clone_sizes.groupby(clone_col).apply(
                        lambda x: (x[x['pre_post'] == 1]['Count'].sum() + pseudo_count) /
                                  (x[x['pre_post'] == 0]['Count'].sum() + pseudo_count))

        after_estimate = clone_sizes.groupby(clone_col).apply(
                        lambda x: (x[x['pre_post'] == 1]['Count'].sum() + pseudo_count))
        before_estimate = clone_sizes.groupby(clone_col).apply(
                        lambda x: (x[x['pre_post'] == 0]['Count'].sum() + pseudo_count))

        return growth_estimate, clone_sizes, after_estimate, before_estimate


    @staticmethod
    def compare_befaft_kl(clone_sizes, pseudo_val=0.0001):
        """
        Entropy

        :param cells_meta: Information on cells
            pd.DataFrame
                Index:
                    RangeIndex
                Columns:
                    * Name: clone, dtype: object
                    * Name: pre_post, dtype: bool
        :return: KL divergence between the pre-growth
        probability distribution and the post-growth distribution.
        :rtype: float
        """

        pre = clone_sizes[clone_sizes['pre_post'] == 0].set_index('clone')['Count']
        post = clone_sizes[clone_sizes['pre_post'] == 1].set_index('clone')['Count']
        pre_prob = pre/pre.sum()
        post_prob = post/post.sum()

        # Pseudo_val added and entropy will renormalize
        prob = pd.concat((pre_prob, post_prob), axis=1).fillna(pseudo_val)

        return entropy(prob.iloc[0], prob.iloc[1])


    @staticmethod
    def compare_kl_across_samples(kl_vals, samples_meta, f_save):
        """ Will plot the distribution between the two phenotype groups
        with Mann-Whitney test.

        :return:
        """

        df = pd.DataFrame({'KL':kl_vals,
                           'Samples':samples_meta})

        g = sns.FacetGrid(df, hue="Samples", palette="Set1")
        g = (g.map(sns.distplot, "KL", hist=False, rug=True))
        plt.savefig(f_save)
        return

    @staticmethod
    def compare_befaft_enrichment():
        return


    @staticmethod
    def compare_phylogeny_triplets(true_lineage, pred_lineage):
        return


    @staticmethod
    def compare_bef_aft(self):
        """Creates a df that contains information on the number of cells from
        each clone before as well as after.
        :return: df.at[ind, "Dominant Before"] = (full_sim.clone_cell == 1).sum() df.at[ind, "Dominant After"]
        = (full_sim.subsample_new_clone_cell == 1).sum()
        """
        return

    def cluster_compare_before_after(self):
        """Compares the performance of clustering on grouping the same clones
        together. :return:
        """
        return


    @staticmethod
    def plot_cluster(cell_af, cell_meta=None, mt_meta=None, f_save=None):
        """
        Args:
            cell_af:
            cell_meta:
            mt_meta:
            f_save:
        """
        ch.plot_cluster(cell_af, row_meta=cell_meta, col_meta=mt_meta,
                        fsave=f_save, to_col_clust=False, to_z=True)



def main():
    return


if "__name__" == "__main__":
    main()

