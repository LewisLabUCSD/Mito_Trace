#import pymc3 as pm
import numpy as np
from numpy import random
import os
import pandas as pd
from tqdm import tqdm
#from src.config import ROOT_DIR
from sklearn.metrics import roc_curve, average_precision_score
from scipy import interp
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
    # TODO Sweep across het, growth rate, coverage, error, cluster size,
    def __init__(self, default_params_f, sweep_params_f):
        self.sweep_params_f = sweep_params_f
        self.default_params_f = default_params_f
        params = read_config_file(default_params_f)
        #self.params =
        sweep_params = read_config_file(sweep_params_f)
        self.sweep_params = sweep_params
        self.default_params = params

        # Create the yaml files in the directory indicated by local_outdir and prefix in sweep_params_f
        self.outdir = os.path.join(sweep_params["local_outdir"], sweep_params["prefix"])
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)


        # Create a grid.
        # Loop through and set each parameter in the the params dictionary.
        self.params_df = pd.DataFrame(list(ParameterGrid(sweep_params['grid'])))
        params_dict = dict()
        # Save the file with a count number. Associate the number to the df
        for ind, val in self.params_df.iterrows():
            f_name = os.path.join(self.outdir, str(ind)+'.yaml' )
            for name, v in val.iteritems():
                if name == 'dominant_clone_sizes':
                    params["initialize"]["clone_sizes"][1] = v
                    # Make the non-clone sum to 1
                    params["initialize"]["clone_sizes"][0] = 1-sum(params["initialize"]["clone_sizes"][1:])
                    print(params["initialize"]["clone_sizes"])
                    assert(np.abs(sum(params["initialize"]["clone_sizes"]) - 1)<0.0001)
                elif name == 'dominant_het':
                    params["het"][0] = v
                elif name == 'dominant_growth':
                    params["growth"]["binomial"]["rates"][1] = v
                else:
                    params = replace_item(params, name, v)
            write_config_file(f_name, params)
            params_dict[f_name] = params
        return

    def run_sweep(self, subset=None):
        """
        :param subset: types:{int, None} default=None, all files used.
                        if list, list of files to run on, if int, number
                        of files to randomly choose.
        :return:
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


    def plot_sensitivity_and_dropout(self, vars=None):
        #f, ax = plt.subplots(1, 3, figsize=(15, 15))

        # coverages = self.params_df['cov_constant'].unique()
        # hets = self.params_df['dominant_het'].unique()
        # err = self.params_df['het_err_rate'].unique()
        # sizes = self.params_df['het_err_rate'].unique()

        metrics = self.params_df.copy()
        metrics['Sensitivity'] = -1
        metrics['% dropout'] = -1

        # Add these results to self.results, which has the meta information too
        for ind, val in self.params_df.iterrows():
            full_sim = self.sweep_results[ind]
            dropout = full_sim.dropout
            prec_scores = full_sim.prec_scores
            rocs = full_sim.rocs
            metrics.at[ind, 'Sensitivity'] = np.mean(prec_scores)
            metrics.at[ind, '% dropout'] = np.mean(dropout)

        # colors, _ = get_colors("categorical", names=coverages,
        #                        n_colors=len(coverages))
        # Seaborn Factorplot
        g = sns.FacetGrid(data=metrics, col="het_err_rate", hue="cov_constant")
        g.map_dataframe(sns.scatterplot, x="dominant_het", y="Sensitivity")
        g.add_legend()
        g.savefig(os.path.join(self.outdir,'sensitivity.png'))


        g = sns.FacetGrid(data=metrics, col="het_err_rate", hue="cov_constant")
        g.map_dataframe(sns.scatterplot, x="dominant_het", y="% dropout")
        g.add_legend()
        g.savefig(os.path.join(self.outdir, 'dropout.png'))
        return

    def plot_ppv(self):
        return

    def cluster_before_after(self):
        for f in self.sweep_results:
            self.sweep_results[f].cluster()
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



# I can make each variable a class?
# Does this ruin running the MCMC? I don't think so, b/c that format is going to be put in after anyway
class FullSimulation:
    def __init__(self, params_f):
        params = read_config_file(params_f)
        self.n_iter = params['num_iterations']
        self.num_cells = params['num_cells']
        self.params = params
        return
        #for i in self.n_iter:

    def run(self):
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
        Will average metrics over simulations.
        :param group: {'init', 'growth', 'both'} This will indicate to group by
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
        cluster_results = []
        print('clustering')
        for s in tqdm(self.sim.values):
            cluster_results.append(s.cluster(s.combined_cell_af))
            print(len(cluster_results[-1]))
        self.cluster_results = cluster_results
        return




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


class Simulation:
    def __init__(self, params_f):
        if isinstance(params_f, str):
            params = read_config_file(params_f)
        else:
            params = params_f

        self.params = params
        check_required(params, ['initialize', 'num_cells', 'num_mt_positions', 'prefix'])
        self.prefix = params['prefix']
        self.num_mt_positions = params['num_mt_positions']
        self.num_cells = params['num_cells']
        if not os.path.exists(params['local_outdir']):
            os.mkdir(params['local_outdir'])
    #should be external method
    def initialize(self):
        self.init_clone_dict()
        self.init_cell_coverage()
        self.init_cell_af()
        #self.init_clone_mt()

    #should be external method
    def grow(self):
        p = self.params
        type = p["growth"]["type"]
        if  type == "poisson":
            self.grow_poisson(p['growth']['poisson'])
        elif type == "binomial":
            self.grow_binomial(p['growth']['binomial'])
        return

    # Static Method
    @staticmethod
    def clone_counts_to_cell_series(clone_counts):
        clone_counts = np.array(clone_counts)
        num_cells = clone_counts.sum()
        clone_cell = -1 * np.ones(shape=[num_cells, ])


        clone_cell[:clone_counts[0]] = 0
        for ind, val in enumerate(clone_counts[1:]):
            start = clone_counts[:ind + 1].sum()
            end = clone_counts[:ind + 1].sum() + val
            # starts at sum(clone_counts[i-1]) ends at clone_counts[i].sum()
            clone_cell[start:end] = ind + 1

        clone_cell = pd.Series(clone_cell, dtype=int)
        return clone_cell

    def init_clone_dict(self):
        ### Add in potential to overwrite the values

        # Gets the clone dictionary. Should also have clone to mt dict.
        clones = self.params['initialize']['clone_sizes']
        num_cells = self.num_cells

        # Option 1: List of fraction of size of each clone. 0s are nonclone size, listed first
        if type(clones) == list:
            #clone_cell = pd.Series(index=range(num_cells))
            clone_counts = np.random.multinomial(num_cells, clones)
            clone_cell  = self.clone_counts_to_cell_series(clone_counts)
            self.clone_cell = clone_cell
        # Option 2: 1 clone. ID'd as 1
        elif type(clones) == int: #One number for dominant clone. the others are not.
            clone_cell = np.zeros(shape=[num_cells,])
            clone_cell[:num_cells] = 1
            clone_cell = clone_cell[::-1]
            clone_cell =  pd.Series(clone_cell, dtype=int)
            self.clone_cell = clone_cell

        # Option 3 To ADD, beta binomial and more complex distributions

        self.num_clones =  len(set(clone_cell.values))-1 # Remove the non-clone
        return clone_cell


    def init_cell_coverage(self):
        """
        There are different modes to the coverage, either a constant or through a distribution.
        :return:
        """
        p = self.params['initialize']['coverage']
        type = p['type']

        num_cells = self.num_cells
        num_pos = self.num_mt_positions
        c = np.zeros([num_cells, num_pos])

        if type == 'constant':
            c[:, :] = p['cov_constant']
        elif type == "poisson":
            # Get the number of coverage per cell based on poisson (should be reads)
            mu_cov_per_cell = p['mu_cov_per_cell']
            num_reads_per_cell = random.poisson(lam=mu_cov_per_cell,
                                                size=num_cells)

            # Number of reads at each position, based on the average for each cell
            for i in num_cells:
                c[i, :] = random.poisson(num_reads_per_cell[i],
                                         size=num_pos)
        self.cells_mt_coverage = c
        return c


    def init_cell_af(self):
        """
        Initialize the cell-by-mtPos af dataframe. Unless a clone:mt dict was provided,
        the first N MT positions will be the clone AFs.
        Creates self.clone_mt_dict and self.cell_af"""

        p = self.params['initialize']

        hets = self.params['het']
        q = self.params['het_err_rate']
        clone_df = self.clone_cell
        num_clones = self.num_clones
        n_cells = self.num_cells
        n_mt = self.num_mt_positions

        # Output
        cell_af = pd.DataFrame(np.zeros(shape=[n_cells, n_mt]))


        if 'mt_clone_map' in p and p['mt_clone_map'] is not None:
            self.clone_mt_dict = p['mt_clone_map']
        else:
            # Each clone points to a mt position
            self.clone_mt_dict = dict()
            for i in range(1,num_clones+1):
                self.clone_mt_dict[i] = i

        # TODO Add the MT clone map so it can contain multiple mutants in lineages

        # If there is a heteroplasmy table in params, it is list of mutant heteroplasmy AFs.
        # If not, will randomly draw based on number of clones
        if type(hets) == list:
            if (len(hets) != num_clones):
                print('here')
            assert(len(hets) == num_clones)

            ## Loop through each clone,
            ## Generate the AF for the clone and non-clones using coverage for each cell
            ## Fill in cell_by_af for that position.
            for ind in range(1, num_clones+1):
                # Generate AF: (clone_df ==  ind).sum()
                n_dom_cells = (clone_df==ind).sum()
                het = hets[ind-1]

                curr_mt = self.clone_mt_dict[ind]


                if p['coverage']['type'] == 'constant':
                    c = p['coverage']['cov_constant']

                    af_i = random.binomial(c, het,
                                           n_dom_cells) / c
                    af_j = random.binomial(c, q,
                                           n_cells - n_dom_cells) / c

                    # Update the dom_cells and non_dom for the current MT
                    cell_af.loc[np.flatnonzero(clone_df == ind), curr_mt] = af_i
                    cell_af.loc[np.flatnonzero(clone_df != ind), curr_mt] = af_j

                # Each cell and position has it's own coverage value, so need to update each
                else:
                    c = self.cells_mt_coverage
                    #Get the cells coverage for the mt position
                    curr_mt_cov= c[:, curr_mt]

                    # Get cell indicies for the clones and nonclones
                    curr_clone_inds = np.flatnonzero(clone_df==ind)
                    curr_nonclone_inds = np.flatnonzero(clone_df!=ind)
                    for cell in curr_clone_inds:
                        # Get one value for curr_mt and cell based on coverage
                        cell_af.loc[cell, curr_mt] = random.binomial(curr_mt_cov[cell], het)
                    for cell in curr_nonclone_inds:
                        cell_af.loc[cell, curr_mt] = random.binomial(curr_mt_cov[cell], q)
                # Loop through each coverage
                #for c in n_dom_cells:

        #####
        # TODO
        # Add noise to the other non-lineage positions
        #####
        self.cell_af = cell_af
        return


    def init_clone_mt(self):
        p = self.params
        if p["initialize"]['type'] == 'growth':
            ## TODO
            # Create a phylogeny and then get the averages of the mutants
            self.average_clone_mt()
        # If not growth, should aready be there.
        return

    def average_clone_mt(self):
        return

    def extract_clone_cells(self, clone_id):
        ids = np.flatnonzero(self.clone_cell == clone_id)
        return ids

    def simulate_expand_cells_af(self, af, growth_inds, sigma):
        """
        Given a cell-by-af vector, expand the AF
        :param af:
        :param growth:
        :param sigma:
        :return:
        """

        new_af = af.iloc[growth_inds].copy() + random.normal(0, sigma, size=af.iloc[growth_inds].shape)
        new_af.index = np.arange(af.index[-1]+1, af.index[-1]+1+new_af.shape[0])
        new_af = pd.concat((af,new_af), axis=0)
        #new_af = np.append(af, np.concatenate(new_af))

        return new_af

    def grow_binomial(self, p):
        timesteps = p["time_steps"]
        rates = p["rates"]

        sigma = self.params['growth']["mutant_af_sigma_noise"]
        cell_af = self.cell_af
        clone_mt_dict = self.clone_mt_dict

        num_clones = self.num_clones+1
        new_dict = {}
        for curr_clone in range(num_clones):
            curr_rate = rates[curr_clone]
            ids = self.extract_clone_cells(curr_clone)
            new_cells = cell_af.loc[ids].copy()
            for i in range(timesteps):
                # Simulate growth for each clone separately.
                growth_inds = np.flatnonzero(random.binomial(1, curr_rate, size=new_cells.shape[0]))
                #new_ids =
                new_cells = self.simulate_expand_cells_af(new_cells, growth_inds, sigma)

            new_dict[curr_clone] = new_cells
            # Create list of cells

        ####TODO
        ## new_lineage_mutants chances. This will see if a mutation will change


        ####TODO
        ## Add death + stimulation rate as well as growth
        # Save the new cell clones df and cell af
        clone_counts = [i.shape[0] for i in new_dict.values()]
        self.new_clone_cell = self.clone_counts_to_cell_series(clone_counts)

        self.new_cell_af = pd.DataFrame(new_dict[0])
        for clone in range(1, self.num_clones+1):
            self.new_cell_af = pd.concat((self.new_cell_af, new_dict[clone]),axis=0).reset_index(drop=True)
        return


    def grow_poisson(self):
        # TODO growth of poisson refactor
        return


    def subsample_new(self, to_delete=False):
        new_cell_af = self.new_cell_af
        p = self.params
        if 'sequence_subsample' in p and p['sequence_subsample'] is not None:
            self.subsample_new_cell_af = new_cell_af.sample(n=self.params['sequence_subsample'])
        else:
            self.subsample_new_cell_af = new_cell_af.sample(n=self.num_cells)

        self.subsample_new_clone_cell = self.new_clone_cell.loc[
            self.subsample_new_cell_af.index]

        if to_delete:
            self.new_cell_af = None
            self.new_clone_cell = None


    def combine_init_growth(self):
        clones = pd.concat(
            (self.clone_cell, self.subsample_new_clone_cell)).reset_index(
            drop=True)
        combined_cell_af = self.cell_af.append(self.subsample_new_cell_af).reset_index(drop=True)

        combined_meta = np.concatenate((np.ones(shape=[self.cell_af.shape[0],]), np.zeros(shape=[self.subsample_new_cell_af.shape[0]])))
        combined_meta = pd.Series(combined_meta, name='After Growth', dtype=int)
        assert(combined_meta.shape[0] == self.cell_af.shape[0]+self.subsample_new_cell_af.shape[0])
        assert (combined_cell_af.shape[0] == self.cell_af.shape[0] +
                self.subsample_new_cell_af.shape[0])
        assert(combined_meta.shape[0] == clones.shape[0])
        assert (combined_cell_af.shape[0] == clones.shape[0])
        self.combined_meta = combined_meta
        self.combined_clones = clones
        self.combined_cell_af = combined_cell_af
        return

    def save(self, f_save=None):
        if f_save is None:
            f_save = os.path.join(self.params['local_outdir'], self.params['prefix']+'.p')
        f = open(f_save, 'wb')
        pickle.dump(self.__dict__, f, 2)
        f.close()

    def save_to_mgatk_format(self):
        """
        Converts into the proper files needed for mgatk. (i.e variant and coverage files)
        :return:
        """

    def load(self):
        filename = self.params['filename']
        f = open(filename, 'rb')
        tmp_dict = pickle.load(f)
        f.close()
        self.__dict__.update(tmp_dict)

    def compare_before_after(self):
        return

    @staticmethod
    def plot_cluster(cell_af, cell_meta=None, mt_meta=None, f_save=None):
        ch.plot_cluster(cell_af, row_meta=cell_meta, col_meta=mt_meta,
                        fsave=f_save, to_col_clust=False, to_z=True)

    @staticmethod
    def cluster(cell_af):
        """
        Dynamic tree clustering of the rows of cell_af
        :param cell_af:
        :return:
        """
        distances = pdist(cell_af, "euclidean")
        link = linkage(distances, "average")
        clusters = cutreeHybrid(link, distances)['labels']
        return clusters

    @staticmethod
    def cluster_kmeans(cell_af):
        distortions = []
        inertias = []
        mapping1 = {}
        mapping2 = {}
        K = range(1, 10)

        for k in K:
            # Building and fitting the model
            kmeanModel = KMeans(n_clusters=k).fit(cell_af)
            kmeanModel.fit(cell_af)

            distortions.append(sum(
                np.min(cdist(cell_af, kmeanModel.cluster_centers_, 'euclidean'),
                       axis=1)) / cell_af.shape[0])
            inertias.append(kmeanModel.inertia_)

            mapping1[k] = sum(
                np.min(cdist(cell_af, kmeanModel.cluster_centers_, 'euclidean'),
                       axis=1)) / cell_af.shape[0]
            mapping2[k] = kmeanModel.inertia_


def main():
    return


if "__name__" == "__main__":
    main()