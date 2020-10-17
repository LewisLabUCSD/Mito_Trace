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

