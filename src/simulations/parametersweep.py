#import pymc3 as pm
import numpy as np
import os
import pandas as pd
import pickle
import seaborn as sns
from pandarallel import pandarallel
import matplotlib.pyplot as plt

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
        params = read_config_file(default_params_f)
        if 'n_clust' not in params:
            self.params = None
        sweep_params = read_config_file(sweep_params_f)
        self.sweep_params_f = sweep_params_f
        self.default_params_f = default_params_f
        self.sweep_params = sweep_params
        self.default_params = params
        self.save_sim = sweep_params['save_sim']


        print(self.save_sim)
        # Create the yaml files in the directory indicated by local_outdir and prefix in sweep_params_f
        self.outdir = os.path.join(sweep_params["outdir"], sweep_params["prefix"])
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        self.f_save = \
            os.path.join(self.outdir,self.sweep_params['prefix'] + '.p')

        #self.tmp_f_save = self.f_save.replace('.p', '') + '_tmp.p'

        self.data_outdir = os.path.join(sweep_params['data_outdir'], sweep_params["prefix"])
        if not os.path.exists(self.data_outdir):
            os.makedirs(self.data_outdir)

        # Create a grid.
        # Loop through and set each parameter in the the params grid parameter.
        self.params_df = pd.DataFrame(list(ParameterGrid(sweep_params['grid'])))
        params_dict = dict()

        for ind, val in self.params_df.iterrows():
            f_name = os.path.join(self.outdir, str(ind)+'.yaml' )
            # Create the name
            params['name'] = str(ind)
            params['data_outdir'] = self.data_outdir
            params['local_outdir'] = self.outdir
            params['prefix'] = self.sweep_params["prefix"]
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

    def patient_samples_simulate(self):
        return

    #dominant_clone_sizes

    def patient_samples_compare(self):
        return

    @staticmethod
    def run_single_sweep(f, outdir, save=True):
        print(f'Running file {f}')
        params_f = os.path.join(outdir, str(f) + '.yaml')
        sim = FullSimulation(params_f)
        sim.run()
        sim.run_metrics()
        if save:
            sim.save()
        return sim.metrics

    def run_sweep(self, subset=None):
        """ Loops through params_df and runs the simulation on that parameter.

        Each run will be saved, and not stored. However, the metrics
        will be stored here.
        :param subset:  Which files to run the simulation. If list,
        list of files to run on, if int, number of files to randomly
        choose. If None, run on all.
        :type subset: int or None or list (default=None)
        """
        params_df = self.params_df
        if isinstance(subset, int):
            params_df = params_df.sample(n=subset)

        sweep_results_df = pd.Series(index=params_df.index, data=params_df.index)
        ###
        #sweep_results_df = sweep_results_df.apply(self.run_single_sweep, args=(self.outdir,))
        pandarallel.initialize(nb_workers=self.sweep_params['cpus'])
        sweep_results_df = sweep_results_df.parallel_apply(self.run_single_sweep, args=(self.outdir, self.save_sim))
        ###

        self.sweep_results = sweep_results_df
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
            curr_results = self.sweep_results[ind]
            dropout = curr_results['dropout']
            prec_scores = curr_results['prec_scores']
            rocs = curr_results['rocs']

            metrics.at[ind, 'Avg. Precision'] = np.mean(prec_scores)
            metrics.at[ind, '% dropout'] = np.mean(dropout)

        metrics = metrics.rename({'cov_constant': 'coverage',
                        'dominant_clone_sizes': 'CHIP proportion',
                        'het_err_rate': 'Error rate',
                        'dominant_het': 'Heteroplasmy'},
                       axis=1)
        self.metrics = metrics

        """df that contains simulation performance metric.
        """
        # Seaborn Factorplot
        g = sns.FacetGrid(data=metrics, col='Error rate',
                          row='CHIP proportion',
                          hue="coverage")
        g.map(sns.stripplot, "Heteroplasmy", "Avg. Precision")
        g.add_legend()
        g.savefig(os.path.join(self.outdir,'precision.png'))

        g = sns.FacetGrid(data=metrics, col='Error rate',
                          row='CHIP proportion',
                          hue="coverage")
        g.map(sns.stripplot, "Heteroplasmy", "% dropout")
        g.add_legend()
        g.savefig(os.path.join(self.outdir, 'dropout.png'))
        return None

    def plot_ppv(self):
        return

    def cluster_before_after(self):
        for f in self.sweep_results:
            self.sweep_results[f].cluster_before_after()
        return

    def plot_before_after_all_clones(self):
        """
        Plot the growth of each clone. Boxplots with x-axis being
        the growth rate for the dominant clone, the y-axis the
        after/before ratio, the color the dominant clone cell size.

        # Expand params_df to include each iteration, which will
        # allow for violinplots. This will lead to
        params_df.shape[0]*num_iterations rows.
        :return:
        """
        df_full = []
        for ind, val in self.params_df.iterrows():
            curr_results = self.sweep_results[ind]
            b_a_df = curr_results['b_a_df']
            # Each element will be a list of the values.
            s = self.params_df.loc[ind].copy()
            curr = pd.DataFrame([s.copy()] * len(b_a_df))
            curr["A/B"] = b_a_df["A/B"].values
            df_full.append(curr)

        df_full = pd.concat(df_full)
        df_full = df_full.astype(
            {'dominant_growth': str, 'dominant_clone_sizes': str,
             'A/B': float})

        sns.violinplot(data=df_full, y="A/B", hue="dominant_growth",
                       x="dominant_clone_sizes")
        plt.savefig(os.path.join(self.outdir, 'growth_before_after.png'))
        #g.savefig(os.path.join(self.outdir, 'precision.png'))
        return None

    def plot_before_after_true_growth(self):
        """
        Plot the growth of cells in each clone before and after.
        Dotplot/heatmap, where each row is a het value, each column
        is growth rate, and the value is dominant clone cells
        growth (after/before). Can make multiple where each column is
        coverage, each row is error rate.

        :return:
        """

        df_full = []
        for ind, val in self.params_df.iterrows():
            curr_results = self.sweep_results[ind]
            b_a_df = curr_results['b_a_df']

            # Each element will be a list of the values.
            s = self.params_df.loc[ind].copy()
            curr = pd.DataFrame([s.copy()] * len(b_a_df))
            curr["A/B"] = b_a_df["A/B"].values
            curr["True clones"] = 1

            curr_pred = curr.copy()
            curr_pred["A/B"] = curr_results['pred_growth_rates']
            curr_pred["True clones"] = 0
            #curr['Predicted A/B'] = curr_results['pred_growth_rates']
            df_full.append(curr.append(curr_pred))

        df_full = pd.concat(df_full)
        df_full = df_full.astype(
            {'dominant_growth': str, 'dominant_clone_sizes': str,
             'A/B': float})

        sns.violinplot(data=df_full, y="A/B", hue="True clones",
                       x="dominant_clone_sizes", row='dominant_growth')
        plt.savefig(os.path.join(self.outdir, 'growth_before_after_pred.png'))
        #g.savefig(os.path.join(self.outdir, 'precision.png'))
        return

    def plot_before_after_cluster_growth(self):
        """
        Plot the growth of cells in each clone before and after.
        Dotplot/heatmap, where each row is a het value, each column
        is growth rate, and the value is dominant clone cell cluster
        growth (after/before). Can make multiple where each column is
        coverage, each row is error rate.
        :return:
        """
        return


    def plot_error_cluster_growth(self):
        """
        Plot the deviation of the growth rates for each pipeline run.
        If A=true growth rate, B=predicted, (C = A-B)
        :return:
        """

    def delete_tmp(self):
        if os.path.exists(self.tmp_f_save):
            os.remove(self.tmp_f_save)

    def save(self, f_save=None):
        if f_save is None:
            f_save = self.f_save
        f = open(f_save, 'wb')
        pickle.dump(self.__dict__, f, 2)
        f.close()

    def load(self, f_save=None):
        if f_save is None:
            f_save = self.f_save
        #filename = self.params['filename']
        f = open(f_save, 'rb')
        tmp_dict = pickle.load(f)
        f.close()
        self.__dict__.update(tmp_dict)