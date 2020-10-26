#import pymc3 as pm
import numpy as np
from numpy import random
import os
import pandas as pd
from tqdm import tqdm
#from src.config import ROOT_DIR
from sklearn.metrics import roc_curve, average_precision_score
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import glob
from sklearn.cluster import KMeans
from sklearn import metrics
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

from .simulation import Simulation

""" Run the simulation similar to in Extended Data Fig 3 from 
    Massively parallel single-cell mitochondrial DNA genotyping and chromatin profiling"""


# I can make each variable a class?
# Does this ruin running the MCMC? I don't think so, b/c that format is going to be put in after anyway
class FullSimulation:
    """
    Class that simulates cell growth for lineage tracing. Reads in a
    parameter file and runs a certain number of iterations based on
    the num_iterations parameter.

    :ivar n_iter: Number of iterations
    :type n_iter: int
    :ivar num_cells: Number of cells to sequence
    :type num_cells: int

    :ivar sim: Each index is a different iteration of the simulation.
    :type sim: pandas Series

    """
    def __init__(self, params_f):
        params = read_config_file(params_f)
        self.n_iter = params['num_iterations']
        self.num_cells = params['num_cells']
        self.params = params
        return
        #for i in self.n_iter:

    def run(self):
        """
        Runs the simulation and stores it in sim attr.

        This uses Pandaralel to parallelize the runs.
        :return:
        """
        # Parallelize df
        df = pd.Series(index=range(self.n_iter))
        df = df.apply(self.run_sim, args=(self.params,))
        #df = df.parallel_apply(self.run_sim, args=(self.params,))

        self.sim = df
        #self.cluster_before_after()
        self.sim_performance_dominant(group='both')
        return

    @staticmethod
    def run_sim(x, params):
        """
        For a simulation, it will initialize, grow, subsample,
        and merge the before stimulus and after stimulus variables.
        :param x: Placeholder variable
        :param params: The parameter dictionary to use
        :return:
        """
        s = Simulation(params)
        s.initialize()
        s.grow()
        s.subsample_new(to_delete=True)
        s.combine_init_growth()
        return s

    def flatten_sim(self):
        ## TODO
        # This will extract out the classes of df
        return

    def sim_performance_dominant(self, group='both'):
        """
        Will colect metrics that are averaged over the simulations.
        These are specifically for looking at the main, dominant clone,
        and what the allele-frequency of that clone variant
        is for each cell.

        :param group: {'init', 'growth', 'both'} This will indicate to group by
        :ivar dropout: Number of dominant clone cells that have 0 reads
        at the lineage variant position.
        :type dropout: list
        :ivar prec_scores: sklearn average precision score based on
        the allele frequencies seen in the dominant clone cells versus
        the non-clone cells.
        :type prec_scores: list
        :ivar rocs: ROC curves for each iteration based on allele
        frequencies.
        :return:
        """
        dropout = []
        rocs = []
        prec_scores = []

        for iter, s in enumerate(self.sim.values):
            # First get the dominant clone , which is indexed as 1
            mt_pos = s.clone_mt_dict[1]
            # TODO account for mt_pos being a list not an int
            if group == 'init':
                clones = s.clone_cell
                cell_af = s.cell_af.loc[:,mt_pos]
            elif group == 'growth':
                clones = s.new_clone_cell
                cell_af = s.new_cell_af.loc[:,mt_pos]
            elif group == 'both':
                #clones = pd.concat((s.clone_cell, s.subsample_new_clone_cell)).reset_index(drop=True)
                #cell_af = pd.concat((s.cell_af.loc[:,mt_pos], s.subsample_new_cell_af.loc[:,mt_pos])).reset_index(drop=True)
                clones = s.combined_clones
                cell_af = s.combined_cell_af.loc[:,mt_pos]
            else:
                raise ValueError('group variable not properly set.')

            y_true = clones.values.copy()
            y_true[y_true != 1] = 0  # Set nondominant clones to 0
            rocs.append(roc_curve(y_true, cell_af))
            prec_scores.append(average_precision_score(y_true, cell_af))
            dropout.append((cell_af[clones==1]==0).sum()/cell_af.shape[0])

        self.dropout = dropout
        self.prec_scores = prec_scores
        self.rocs = rocs
        return


    def reduce_cells(self, cell_af):
        #self.sim
        return


    def cluster_before_after(self):
        """
        Loops through the simulations and for each,
        it clusters the cells.

        :ivar cluster_results: Cluster labels for each cell in each
        iteration.
        :type List of tuples, which is a list of
        a tuple, where the tuple is indexed by the cell and the value
        is the cell's cluster label
        """
        cluster_results = []
        print('clustering')
        for s in tqdm(self.sim.values):
            cluster_results.append(s.cluster(s.combined_cell_af))
            print(len(cluster_results[-1]))
        self.cluster_results = cluster_results

    def save(self, f_save=None):
        if f_save is None:
            f_save = os.path.join(self.params['local_outdir'], self.params['prefix']+'.p')
        f = open(f_save, 'wb')
        pickle.dump(self.__dict__, f, 2)
        f.close()

    def load(self, filename):
        #filename = self.params['filename']
        f = open(filename, 'rb')
        tmp_dict = pickle.load(f)
        f.close()
        self.__dict__.update(tmp_dict)


def main():
    return


if "__name__" == "__main__":
    main()