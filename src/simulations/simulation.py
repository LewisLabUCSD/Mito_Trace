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


class Simulation:
    """Lineage tracing simulation of one sample

    Will initialize cells based on their parameters and grow as well. This
    should be a flexible framework, to add different ways to initialize, grow,
    and metrics to have. Additionally can cluster these results.

    :ivar params
    :type params: dict
    """

    def __init__(self, params_f):
        """
        :param params_f: Parameter yaml file for the specifications
        :type params_f: yaml file or dict
        """
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
        """
        Args:
            clone_counts:
        """
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
        """There are different modes to the coverage, either a constant or
        through a distribution. :return:
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
        """Initialize the cell-by-mtPos af dataframe. Unless a clone:mt dict was
        provided, the first N MT positions will be the clone AFs. Creates
        self.clone_mt_dict and self.cell_af
        """

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
        """
        Args:
            clone_id:
        """
        ids = np.flatnonzero(self.clone_cell == clone_id)
        return ids

    def simulate_expand_cells_af(self, af, growth_inds, sigma):
        """Given a cell-by-af vector, expand the AF.

        Expanded AF occurs by duplicating cells that grew based on the
        growth_inds vector. It will add standard error to each af based on sigma
        :param af: :param growth: Indices of AF to copy :param sigma: Variance
        to add to AF of child. :return:

        Args:
            af:
            growth_inds:
            sigma:
        """

        new_af = af.iloc[growth_inds].copy() + random.normal(0, sigma, size=af.iloc[growth_inds].shape)
        new_af.index = np.arange(af.index[-1]+1, af.index[-1]+1+new_af.shape[0])
        new_af = pd.concat((af,new_af), axis=0)
        #new_af = np.append(af, np.concatenate(new_af))
        return new_af

    def grow_binomial(self, p):
        """
        Args:
            p:
        """
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
        """
        Args:
            to_delete:
        """
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
        assert(combined_cell_af.shape[0] == clones.shape[0])
        self.combined_meta = combined_meta
        self.combined_clones = clones
        self.combined_cell_af = combined_cell_af
        return

    def save(self, f_save=None):
        """
        Args:
            f_save:
        """
        if f_save is None:
            f_save = os.path.join(self.params['local_outdir'], self.params['prefix']+'.p')
        f = open(f_save, 'wb')
        pickle.dump(self.__dict__, f, 2)
        f.close()

    def save_to_mgatk_format(self):
        """Converts into the proper files needed for mgatk. (i.e variant and
        coverage files) :return:
        """

    def load(self):
        filename = self.params['filename']
        f = open(filename, 'rb')
        tmp_dict = pickle.load(f)
        f.close()
        self.__dict__.update(tmp_dict)

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
    def cluster_kmeans(cell_af):
        """
        Args:
            cell_af:
        """
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