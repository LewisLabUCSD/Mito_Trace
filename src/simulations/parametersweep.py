#import pymc3 as pm
import numpy as np
import os
import pandas as pd
import pickle
import seaborn as sns
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


from .fullsimulation import FullSimulation

""" Run the simulation similar to in Extended Data Fig 3 from 
    Massively parallel single-cell mitochondrial DNA genotyping and chromatin profiling"""


def replace_item(obj, key, replace_value):
    # https: // stackoverflow.com / a / 45335542
    for k, v in obj.items():
        if isinstance(v, dict):
            obj[k] = replace_item(v, key, replace_value)
    if key in obj:
        obj[key] = replace_value
    return obj


class ParameterSweep:
    """
    A class for running a lineage simulation with varying parameters.

    :ivar outdir: Directory to save all output to
    :ivar sweep_params: The hyperparameters dictionary used
    :ivar default_params: The baseline parameters used.
    :ivar metrics: pd DataFrame that contains output metrics
    """
    def __init__(self, default_params_f, sweep_params_f):
        """

        Initialize ParameterSweep creates the pipeline yaml files to be
        run with run_sweep.
        :param default_params_f: File that contains the baseline
        parameters for the simulation that all runs will use.
        Grid is the key in the parameter file that contains the
        hyperparameters.
        :type default_params_f: str
        :param sweep_params_f: File that contains the hyperparameters to
        do a grid of all parameters to run the pipeline on.
        :type str

        """
        self.sweep_params_f = sweep_params_f
        self.default_params_f = default_params_f
        params = read_config_file(default_params_f)
        sweep_params = read_config_file(sweep_params_f)
        self.sweep_params = sweep_params
        self.default_params = params

        # Create the yaml files in the directory indicated by local_outdir and prefix in sweep_params_f
        self.outdir = os.path.join(sweep_params["local_outdir"], sweep_params["prefix"])
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Create a grid.
        # Loop through and set each parameter in the the params grid parameter.
        self.params_df = pd.DataFrame(list(ParameterGrid(sweep_params['grid'])))
        params_dict = dict()

        for ind, val in self.params_df.iterrows():
            f_name = os.path.join(self.outdir, str(ind)+'.yaml' )
            for name, v in val.iteritems():
                # Set the specific variables that need to be updated
                if name == 'dominant_clone_sizes':
                    params["initialize"]["clone_sizes"][1] = v
                    # Set the non-clone such that all clones sum to 1
                    params["initialize"]["clone_sizes"][0] = 1-sum(params["initialize"]["clone_sizes"][1:])
                    assert(np.abs(sum(params["initialize"]["clone_sizes"]) - 1)<0.0001)
                elif name == 'dominant_het':
                    params["het"][0] = v
                elif name == 'dominant_growth':
                    params["growth"]["binomial"]["rates"][1] = v
                else: # These parameters assumed to have the same name
                      # In both files
                    params = replace_item(params, name, v)
            write_config_file(f_name, params)
            params_dict[f_name] = params
        return

    def run_sweep(self, subset=None):
        """
        Loops through self.params_df and runs the simulation on that
        parameter.

        :param subset:  Which files to run the simulation. If list,
        list of files to run on, if int, number of files to randomly
        choose. If None, run on all.
        :type subset: int or None or list (default=None)

        """

        params_df = self.params_df
        if isinstance(subset, int):
            params_df = params_df.sample(n=subset)

        sweep_results = dict()
        for f, val in params_df.iterrows():
            params_f = os.path.join(self.outdir, str(f) +'.yaml')
            print(f"Running with file: {f}")
            sim = FullSimulation(params_f)
            sim.run()
            sweep_results[f] = sim
        self.sweep_results = sweep_results
        return


    def plot_sensitivity_and_dropout(self):
        """
        Scatterplots of the average precision score and dropout
        against variant heteroplasmy.

        Additional simulation parameters used for groupings are
        coverage (color), and  error rate (column).
        Uses sklearn's average_precision_score to get precision.
        For dropout, estimates how many times the counts are 0 in the
        clone mitochondrial variant.
        """

        metrics = self.params_df.copy()
        metrics['Avg. Precision'] = -1.0
        metrics['% dropout'] = -1.0

        # Add these results to self.results, which has the meta information too
        for ind, val in self.params_df.iterrows():
            full_sim = self.sweep_results[ind]
            dropout = full_sim.dropout
            prec_scores = full_sim.prec_scores
            rocs = full_sim.rocs

            metrics.at[ind, 'Avg. Precision'] = np.mean(prec_scores)
            metrics.at[ind, '% dropout'] = np.mean(dropout)

        self.metrics = metrics
        # Seaborn Factorplot
        g = sns.FacetGrid(data=metrics, col="het_err_rate", hue="cov_constant")
        g.map(sns.scatterplot, "dominant_het", "Avg. Precision")
        g.add_legend()
        g.savefig(os.path.join(self.outdir,'precision.png'))

        g = sns.FacetGrid(data=metrics, col="het_err_rate", hue="cov_constant")
        g.map(sns.scatterplot, "dominant_het", "% dropout")
        g.add_legend()
        g.savefig(os.path.join(self.outdir, 'dropout.png'))
        return None


    def plot_ppv(self):
        return

    def cluster_before_after(self):
        for f in self.sweep_results:
            self.sweep_results[f].cluster()
        return

    def plot_before_after_all_clones(self):
        """
        Plot the growth of each clone.
        :return:
        """
        return

    def plot_before_after_true_growth(self):
        """
        Plot the growth of cells in each clone before and after
        :return:
        """
        return

    def plot_before_after_cluster_growth(self):
        """
        Plot the growth of the predicted dominant clone based on clustering
        results and picking the top clone.
        :return:
        """
        return

    def save(self):
        f_save = os.path.join(self.outdir, self.sweep_params['prefix']+'.p')
        f = open(f_save, 'wb')
        pickle.dump(self.__dict__, f, 2)
        f.close()

    def load(self, filename):
        #filename = self.params['filename']
        f = open(filename, 'rb')
        tmp_dict = pickle.load(f)
        f.close()
        self.__dict__.update(tmp_dict)