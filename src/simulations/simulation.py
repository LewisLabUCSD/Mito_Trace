#import pymc3 as pm
import numpy as np
from numpy import random
import os
import pandas as pd
import pickle
from src.simulations.utils.config import read_config_file, write_config_file
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
        """ (1) Pre-growth cell population is instantiated.

        Creates a clone-MT dictionary, cell coverage matrix
        (or an int, depending on parameters), and cell-AF matrix.
        :return:
        """
        self.init_clone_dict()
        self.init_cell_coverage()
        self.init_cell_af()
        #self.init_clone_mt()

    #should be external method
    def grow(self):
        """ (2) Growth of cells is run."""
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
        """ Generates new cell IDs based on cluster count iterable
        :param clone_counts: Each i'th element is the number of cells in
        cluster i.
        :type clone_counts: iterable
        :return each index name is a cell ID and each value is which cluster
        the cell belongs too.
        :rtype pd.Series
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
        """1A
        """

        ### Add in potential to overwrite the values
        # Gets the clone dictionary. Should also have clone to mt dict.
        clones = self.params['initialize']['clone_sizes']

        if 'num_cells_population' not in self.params:
            self.num_cells_pop = self.num_cells
        else:
            self.num_cells_pop = self.params['num_cells_population']

        num_cells = self.num_cells_pop

        # Option 1: List of fraction of size of each clone. 0s are nonclone size, listed first
        if type(clones) == list:
            #clone_cell = pd.Series(index=range(num_cells))
            clone_counts = np.random.multinomial(num_cells, clones)
            clone_cell  = self.clone_counts_to_cell_series(clone_counts)
            self.clone_cell_pop = clone_cell

            # Choose subset to be sampled
            self.clone_cell = clone_cell.sample(n=self.num_cells).sort_values()
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
        """1B

        There are different modes to the coverage, either a constant or
        through a distribution.
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


    @staticmethod
    def create_cell_af(clone_df, mt_dict, n_cells, n_mt, num_clones,
                       cov_params, hets, het_err, coverage=None):
        cell_af = pd.DataFrame(np.zeros(shape=[n_cells, n_mt]))

        #p = self.params['initialize']

        ## Loop through each clone,
        ## Generate the AF for the clone and non-clones using coverage for each cell
        ## Fill in cell_by_af for that position.
        for ind in range(1, num_clones + 1):
            # Generate AF: (clone_df ==  ind).sum()
            n_dom_cells = (clone_df == ind).sum()
            het = hets[ind - 1]

            curr_mt = mt_dict[ind]

            if cov_params['coverage']['type'] == 'constant':
                c = cov_params['coverage']['cov_constant']

                af_i = random.binomial(c, het, n_dom_cells) / c
                af_j = random.binomial(c, het_err, n_cells - n_dom_cells) / c

                # Update the dom_cells and non_dom for the current MT
                cell_af.loc[
                    np.flatnonzero(clone_df == ind), curr_mt] = af_i
                cell_af.loc[
                    np.flatnonzero(clone_df != ind), curr_mt] = af_j

            # Each cell and position has it's own coverage value, so need to update each
            else:
                if coverage is None:
                    raise("coverage needs to be assigned")
                c = coverage

                # Get the cells coverage for the mt position
                curr_mt_cov = c[:, curr_mt]

                # Get cell indicies for the clones and nonclones
                curr_clone_inds = np.flatnonzero(clone_df == ind)
                curr_nonclone_inds = np.flatnonzero(clone_df != ind)
                for cell in curr_clone_inds:
                    # Get one value for curr_mt and cell based on coverage
                    cell_af.loc[cell, curr_mt] = random.binomial(
                        curr_mt_cov[cell], het)
                for cell in curr_nonclone_inds:
                    cell_af.loc[cell, curr_mt] = random.binomial(
                        curr_mt_cov[cell], het_err)
        return cell_af


    ##########
    def init_cell_af(self):
        """1C. Initialize the cell-by-mtPos af dataframe. Unless a clone:mt dict was
        provided, the first N MT positions will be the clone AFs. Creates
        self.clone_mt_dict and self.cell_af
        """
        p = self.params['initialize']
        hets = self.params['het']
        clone_df = self.clone_cell
        num_clones = self.num_clones
        n_cells = self.num_cells
        n_mt = self.num_mt_positions

        # Get the MT map
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
            assert(len(hets) == num_clones)

        # Get the cell_af based on MT dictionary and cell coverage
        self.cell_af = self.create_cell_af(clone_df, self.clone_mt_dict,
                                           n_cells, n_mt, num_clones,
                                           self.params['initialize'],
                                           hets,
                                           self.params['het_err_rate'],
                                           coverage=None)
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

    @staticmethod
    def extract_clone_cells(clone_cell, clone_id):
        """ Returns the numbered indices of the specific clones

        :param clone_cell: Each element is the indexed cell's clone label.
        :type clone_cell: np array or pd.Series

        :param clone_id:
        :type clone_id: int or string
        """
        ids = np.flatnonzero(clone_cell == clone_id)
        return ids

    @staticmethod
    def simulate_expand_cells_af(af, growth_inds, sigma):
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
        """ (2.1.2)
        :param p: contains time_steps, rates,
        :type dict
        """
        timesteps = p["time_steps"]
        rates = p["rates"]

        num_clones = self.num_clones+1
        new_dict = {}
        for curr_clone in range(num_clones):
            curr_rate = rates[curr_clone]
            ids = self.extract_clone_cells(self.clone_cell_pop, curr_clone)
            num_curr_cells = len(ids)

            for i in range(timesteps):
                # Simulate growth for each clone separately.
                growth_inds = (random.binomial(1, curr_rate, size=num_curr_cells)).sum()
                num_curr_cells += growth_inds.sum()

            new_dict[curr_clone] = num_curr_cells


        ####TODO
        ## new_lineage_mutants chances. This will see if a mutation will change
        ####TODO
        ## Add death + stimulation rate as well as growth
        # Save the new cell clones df and cell af
        clone_counts = [i for i in new_dict.values()]
        self.new_clone_cell = self.clone_counts_to_cell_series(clone_counts)
        # Do not make cell_af, will make this only when subsampled.


        # self.new_cell_af = pd.DataFrame()
        # for clone in range(1, self.num_clones+1):
        #     self.new_cell_af = pd.concat((self.new_cell_af, new_dict[clone]),axis=0).reset_index(drop=True)
        return


    def grow_binomial_old(self, p):
        """ (2.1.1)
        :param p: contains time_steps, rates,
                 and [growth][mutant_af_sigma_noise
        :type dict
        """
        timesteps = p["time_steps"]
        rates = p["rates"]

        sigma = self.params['growth']["mutant_af_sigma_noise"]
        cell_af = self.cell_af

        num_clones = self.num_clones+1
        new_dict = {}
        for curr_clone in range(num_clones):
            curr_rate = rates[curr_clone]
            ids = self.extract_clone_cells(self.clone_cell, curr_clone)
            new_cells = cell_af.iloc[ids].copy()
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


    def grow_poisson(self, p):
        # TODO growth of poisson refactor
        return


    def subsample_new(self, to_delete=False):
        """(3) Subsample from new cell population and generate cell_af

        :param to_delete: To remove the cells that grew (which takes up
        a lot of RAM).
        :type to_delete: bool
        """
        p = self.params
        if 'sequence_subsample' in p and p['sequence_subsample'] is not None:
            self.subsample_new_clone_cell = self.new_clone_cell.sample(n=self.params['sequence_subsample'])
        else:
            self.subsample_new_clone_cell = self.new_clone_cell.sample(
                n=self.num_cells)

        #print(f'New cell af, {len(self.subsample_new_clone_cell)}')
        # Generate subsample_new_cell_af
        self.subsample_new_cell_af = self.create_cell_af(clone_df=self.subsample_new_clone_cell,
                                                         mt_dict = self.clone_mt_dict,
                                                         n_cells=len(self.subsample_new_clone_cell),
                                                         n_mt=self.num_mt_positions,
                                                         num_clones=self.num_clones,
                                                         cov_params=p['initialize'],
                                                         hets=
                                                             self.params[
                                                                 'het'],
                                                         het_err=self.params['het_err_rate'],
                                                         coverage=None
                                                         )

        if to_delete:
            self.new_cell_af = None
            self.new_clone_cell = None


    def subsample_new_old(self, to_delete=False):
        """(3) Subsample from new cell population

        :param to_delete: To remove the cells that grew (which takes up
        a lot of RAM).
        :type to_delete: bool
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
        """(4) Add the pre- and post- population of cells into a group.

        :return:
        """
        combined_cell_af = self.cell_af.append(self.subsample_new_cell_af).reset_index(drop=True)
        combined_clones = pd.concat(
            (self.clone_cell, self.subsample_new_clone_cell)).reset_index(
            drop=True)

        combined_befaft = np.concatenate((np.zeros(shape=[self.cell_af.shape[0],]), np.ones(shape=[self.subsample_new_cell_af.shape[0]])))
        combined_meta = pd.DataFrame({"pre_post": combined_befaft, "clone": combined_clones})
        #combined_meta = pd.Series(combined_meta, name='After Growth', dtype=int)
        assert(combined_meta.shape[0] == self.cell_af.shape[0]+self.subsample_new_cell_af.shape[0])
        assert (combined_cell_af.shape[0] == self.cell_af.shape[0] +
                self.subsample_new_cell_af.shape[0])
        assert(combined_meta.shape[0] == combined_clones.shape[0])
        assert(combined_cell_af.shape[0] == combined_clones.shape[0])
        self.combined_meta = combined_meta
        self.combined_clones = combined_clones
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

    @staticmethod
    def expand_to_mgatk(curr_mt_af,mt_ref):
        ref = mt_ref[curr_mt_af.name]
        pos = curr_mt_af.name
        return pd.DataFrame({"Ref":ref, "Pos":pos, "Val":curr_mt_af})

    def test_save_to_mgatk_format(self):
        df = pd.DataFrame( [[10,0,1,3,5], [3,0,5,5,0], [6,2,1,1,0]] , columns=np.arange(0,5))
        mt_ref_dict = {0: "A", 1: "G", 2: "C", 3: "C", 4: "T"}
        mt_ref = pd.DataFrame({"Pos": mt_ref_dict.keys(), "Ref": mt_ref_dict})
        return

    def save_to_mgatk_format(self, mt_ref, out_f):
        """Converts into the proper files needed for mgatk. (i.e variant and
        coverage files)

        :return:
        """
        cell_af = self.subsample_new_cell_af
        chars = ["A", "G", "C", "T"]
        def alt_generate(x):
            curr = chars.copy()
            curr.remove(x["Ref"])
            return np.random.choice(curr)
        alt_ref = mt_ref.apply(alt_generate, axis=1)

        # First use the AF and choose an alternative allele
        df_stack = cell_af.stack().reset_index().rename(
            {"level_0": "Cell", "level_1": "MT_pos", 0: "Coverage"},
            axis=1)
        df_stack["Nucleotide"] = df_stack["MT_pos"].apply(
            lambda x: alt_ref[x])

        # Add on the reference allele
        df_stack_ref = cell_af.stack().reset_index().rename(
            {"level_0": "Cell", "level_1": "MT_pos", 0: "Coverage"},
            axis=1)
        df_stack_ref["Coverage"] = 1-df_stack_ref["Coverage"]
        df_stack["Nucleotide"] = df_stack["MT_pos"].apply(
            lambda x: mt_ref[x])

        # Save the NTs.
        # For concordance, split the coverage in two
        df_stack = pd.concat(df_stack, df_stack_ref)
        for ind, val in df_stack.groupby("Nucleotide"):
            # Drop the 0s
            curr = val[val["Coverage"]>0]
            # Save file
            curr_out_f = out_f + "_" + ind + ".txt"
            curr.to_csv(curr_out_f)

        # Save the coverage.
        coverage = self.cells_mt_coverage
        if type(coverage) != int:
            coverage_stack = pd.DataFrame(coverage).stack().reset_index().rename(
                {"level_0": "Cell", "level_1": "MT Position", 0: "Coverage"},
                axis=1)
        else:
            coverage_stack = pd.DataFrame(self.cells_mt_coverage)*np.ones(shape=cell_af.shape).stack().reset_index().rename(
                {"level_0": "Cell", "level_1": "MT Position",  0: "Coverage"},
                axis=1)
        curr_out_f = out_f + "_" + "coverage.txt"
        coverage_stack.to_csv(curr_out_f)
        return

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


def main():
    return


if "__name__" == "__main__":
    main()