#import pymc3 as pm
import numpy as np
from numpy import random
import os
from tqdm import tqdm
# Use numpy random to generate simulated data
#from src.config import ROOT_DIR
from sklearn.metrics import roc_curve, average_precision_score
from scipy import interp
import matplotlib.pyplot as plt
from mplh.color_utils import get_colors
from mplh.fig_utils import legend_from_color
from src.config import FIGURES_DIR

try:
    from .dominant_clone import *
except ModuleNotFoundError:
    from dominant_clone import *


os.chdir(FIGURES_DIR)
outdir = "simulation/extfig3/"
if not os.path.exists(outdir):
    os.makedirs(outdir)
f_save = os.path.join(outdir, "extFig3_ij")

n_sim = 10000
n_cells = 1000
n_dom_cells = 100
coverages= [10, 20, 50]

hets= [0.01, 0.1, 0.2,0.4,0.8]
err_het=0.1



""" Run the simulation similar to in Extended Data Fig 3 from 
    Massively parallel single-cell mitochondrial DNA genotyping and chromatin profiling"""
def simulate(n_sim, n_cells, n_dom_cells, coverage, het, err_het):
    dropout = np.zeros([n_sim, ])
    cell_af = np.zeros([n_sim, n_cells])
    rocs = []
    prec_scores = []
    y_true = np.append(np.ones(n_dom_cells),np.zeros(n_cells-n_dom_cells) )

    for sim in tqdm(range(n_sim)):
        af_i = random.binomial(coverage,het, n_dom_cells)/coverage
        af_j = random.binomial(coverage, err_het, n_cells-n_dom_cells)/coverage

        dropout[sim] = (af_i == 0).sum()
        cell_af[sim, :n_dom_cells] = af_i
        cell_af[sim, n_dom_cells:] = af_j

        rocs.append(roc_curve(y_true, cell_af[sim]))
        prec_scores.append(average_precision_score(y_true, cell_af[sim]))

    # Construct sensitivity-specificity
    #cell_af[sim, y_true]

    return cell_af, y_true, rocs, prec_scores, dropout


def expand_cell_af(cell_af, growth, sigma):
    return random.normal(cell_af, sigma, size=growth)


def simulate_expand_cells_af(af, growth, sigma):
    """
    Given a cell-by-af vector, expand the AF
    :param af:
    :param growth:
    :param sigma:
    :return:
    """
    new_af = []
    for ind, cell in enumerate(af):
        new_af.append(expand_cell_af(cell, growth[ind], sigma ))

    new_af = np.append(af, np.concatenate(new_af))
    return new_af


def grow_poisson(af, clone_meta, clone_growth, non_clone_growth, num_cells, sigma=0.05):
    # First group by clones, then simulate growth, then sample number of cells

    ## ToDO for multiple clones
    #new_af = dict() # Dictonary for each clone
    # clones = set(clone_meta)
    #for c in clones:
        # curr = af[np.where(clone_meta==c)]
    clones_af = af[clone_meta==1]
    nonclones_af = af[clone_meta==0]

    # Sample from poisson the growth
    grow_clones = random.poisson(lam=clone_growth, size=(clone_meta == 1))
    grow_nonclones = random.poisson(lam=non_clone_growth, size=(clone_meta == 0))


    # Expand the cells by adding and changing vary their AFs by normal distribution
    new_c = simulate_expand_cells_af(clones_af, grow_clones, sigma=sigma)
    new_nonc = simulate_expand_cells_af(nonclones_af, grow_nonclones, sigma=sigma)

    # Subsample from the new cells
    new_meta = np.append(np.ones(shape=[len(new_c),]),
                         np.zeros(shape=[len(new_nonc),]))

    # Subsample num_cells by sampling from the clone IDs.
    # Then just take the first indices for the clone and nonclones, since
    # the ordering is random anyway.
    new_meta = random.choice(new_meta, size=num_cells , replace=False)
    new_af = np.append(new_c[:(new_meta==1).sum()],
                       new_nonc[:(new_meta==0).sum()])

    #new_af[c]
    return new_af, new_meta


def simulate_growth(cell_af, clone_meta, clone_growth, non_clone_growth=0.5, num_cells=None):
    """

    :param cell_af: iteration-by-cells AF
    :param clone_meta: |cells| np array indicating if it is the clone or not
    :param clone_growth: Poisson parameter for the clone cells
    :param non_clone_growth: poisson parameter for non-clone cells
    :param num_cells: Number of cells to sample from after growth. If None, defaults to the current number of cells
    :return:
    """
    n_sim = cell_af.shape[0]
    n_cells_orig = cell_af.shape[1]

    # For each clone, see how they grow with lambda
    if num_cells is None:
        num_cells = n_cells_orig

    cell_af_growth = np.zeros(n_sim, num_cells)
    cell_af_meta = -1*np.ones(n_sim, num_cells)

    for sim in tqdm(range(n_sim)):
        new_af, new_meta = grow_poisson(cell_af[sim], clone_meta, clone_growth, non_clone_growth, num_cells)
        cell_af_growth[sim] = new_af
        cell_af_meta[sim] = new_meta

    return cell_af_growth, cell_af_meta


def generate_coverage(num_cells, num_pos=None, cov_constant=None, mu_cov_per_cell=None, mu_dist_reads=None, type='constant'):
    """
    There are different modes to the coverage, either a constant or through a distribution.
    :return:
    """
    c = np.zeros([num_cells, num_pos])
    if type == 'constant':
        c[:, :] = cov_constant
        return c
    elif type=="poisson":
        # Get the number of coverage per cell based on poisson (should be reads)
        num_reads_per_cell = random.poisson(lam=mu_cov_per_cell, size=num_cells)

        # Number of reads at each position, based on the average for each cell
        for i in num_cells:
            c[i,:] = random.poisson(num_reads_per_cell[i], size=num_pos)
    return c


def plot_multiple_rocs(rocs):
    ## Needed for putting all curves in same space
    base_fpr = np.linspace(0, 1, 101)
    tprs = []
    for roc in rocs:
        fpr, tpr = roc[0], roc[1]
        tpr = interp(base_fpr, fpr, tpr)
        tpr[0] = 0.0
        tprs.append(tpr)

    tprs = np.array(tprs)
    print(tprs)
    print(len(tprs))
    mean_tprs = tprs.mean(axis=0)
    std = tprs.std(axis=0)

    tprs_upper = np.minimum(mean_tprs + std, 1)
    tprs_lower = mean_tprs - std

    plt.plot(base_fpr, mean_tprs, 'b')
    plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='grey',
                     alpha=0.3)

    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.axes().set_aspect('equal', 'datalim')
    plt.show()
    return


def extfig3_g(n_sim, n_cells, n_dom_cells, coverage, het, err_het):
    cell_af, y_true, rocs, prec_scores, dropout = simulate(n_sim, n_cells, n_dom_cells, coverage, het, err_het)
    plot_multiple_rocs(rocs)
    return


def extfig3_ij(n_sim, n_cells, n_dom_cells, coverages, hets, err_het, f_save):
    ## Need to add in the threshold for counting

    ## "For ‘detection’, we required the cell to have at least half of the simulated
    ## heteroplasmy (p/2)"
    ## Do they remove the cell from the score if not enough coverage?


    f,ax = plt.subplots(1,3, figsize=(15,15))
    sensitivity_dict = {}
    dropout_dict = {}
    for coverage in coverages:
        for het in hets:
            cell_af, y_true, rocs, prec_scores, dropout = simulate(n_sim, n_cells,
                                             n_dom_cells, coverage, het,
                                             err_het)

            sensitivity_dict[(coverage, het)] = np.mean(prec_scores)
            dropout_dict[(coverage, het)] = np.mean(dropout)

    colors,_ = get_colors("categorical",names=coverages, n_colors=len(coverages))
    #for c in len(coverages)
    for coverage in coverages:
        # Create scatter for hets
        sens = []
        drop = []
        for het in hets:
            sens.append(sensitivity_dict[(coverage, het)])
            drop.append(dropout_dict[(coverage,het)])
        ax[0].scatter(hets,sens, color=colors[coverage])
        ax[1].scatter(hets,drop, color=colors[coverage])

    ax[0].set_ylabel("Sensitivity")
    ax[0].set_xlabel("% Heteroplasmy")
    ax[1].set_ylabel("% dropout")
    ax[1].set_xlabel("% Heteroplasmy")

    plt.suptitle(f"error het: {err_het}\n Clone cells: {n_dom_cells}\n"
                 f"Total Cells: {n_cells}")
    legend_from_color(colors, curr_ax=ax[1])
    plt.show()
    plt.savefig(f_save.replace(".png","")+".png")
    return


extfig3_ij(n_sim, n_cells, n_dom_cells, coverages, hets, err_het, f_save)

