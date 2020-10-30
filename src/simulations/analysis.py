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
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import cdist
from pandarallel import pandarallel
pandarallel.initialize(nb_workers=32)

from mplh.color_utils import get_colors
from mplh.fig_utils import legend_from_color
from mplh import cluster_help as ch
from src.simulations.utils.config import read_config_file, write_config_file

from dynamicTreeCut import cutreeHybrid
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage
from sklearn.model_selection import ParameterGrid
from src.simulations.utils.config import check_required


class Analysis:
    """Analysis of lineage tracing experiments

    Can look at both supervised (from simulation)- where the clones are
    known, along with unsupervised, which is the experimental side.

    :ivar params
    :type params: dict
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
    def cluster_kmeans(X, min_c=2, max_c=10):
        """
        Args:
            cell_af:
        """

        best_n_clusters = -1
        sil_score_max = -1  # this is the minimum possible score
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
        model = KMeans(n_clusters=best_n_clusters, init='k-means++',
                       max_iter=100, n_init=1)
        labels = model.fit_predict(X)

        return labels

    

    def compare_before_after(self):
        """Creates a df that contains information on the number of cells from
        each clone before as well as after. :return: df.at[ind, "Dominant
        Before"] = (full_sim.clone_cell == 1).sum() df.at[ind, "Dominant After"]
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

