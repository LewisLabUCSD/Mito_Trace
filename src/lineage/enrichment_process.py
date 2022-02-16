from os.path import join, exists, dirname
from glob import glob
import pickle
import mplh.cluster_help as ch
import os
#import vireoSNP
import numpy as np
from scipy import sparse
from scipy.io import mmread
import matplotlib.pyplot as plt
from scipy.stats import hypergeom, fisher_exact
#print(vireoSNP.__version__)
import click
import pandas as pd
import seaborn as sns
#from vireoSNP import Vireo
import matplotlib as mpl
mpl.use('Agg')

from icecream import ic
#np.set_printoptions(formatter={'float': lambda x: format(x, '.5f')})

#n_clone_list = [3, 5, 10, 20, 40]  # [2,3,4,5,6,7]

def plot_volcano(enrich_stats, x="Flt3l fold enrichment",
                 y="Fisher -log10p", hue=None, f_save=None, v=0,
                 size=None, to_close=True, to_log=True, ylim=None,
                 xlim=None):
    ic(enrich_stats.head())
    enrich_stats = enrich_stats.astype("float", errors='ignore')

    f, ax = plt.subplots(figsize=(10, 10), sharey=True, sharex=True)

    if to_log:
        enrich_stats[f"log2 {x}"] = np.log2(enrich_stats[x])
        x = f"log2 {x}"
        v=0

    if hue is None:
        n_clr=1
    else:
        n_clr=len(set(enrich_stats[hue]))
    sns.scatterplot(data=enrich_stats, x=x,
                    y=y, s=100, sizes=(20, 200),
                    palette=sns.color_palette("Set1", n_clr),
                    ax=ax, hue=hue, size=size)
    plt.axvline(x=v)
    #ax.plot([0.5], [0.5], transform=ax.transAxes, color='black')
    #plt.axis('square')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
               borderaxespad=0)
    if f_save is not None:
        plt.savefig(f_save, bbox_inches='tight')
    if to_close:
        plt.close()
    return


def wrap_plot_volcano(enrich_stats_all, x="Flt3l fold enrichment norm",
                 y="Fisher -log10p", hue=None, f_save=None, v=0,
                 size=None, to_close=False, to_log=False):

    f, ax = plt.subplots(figsize=(15, 15), ncols=len(enrich_stats_all),
                         sharey=True, sharex=True)

    count = 0
    for k in enrich_stats_all:
        curr_ax = ax[count]

        enrich_stats = enrich_stats_all[k].copy()
        enrich_stats = enrich_stats.astype(float)

        if to_log:
            new_x = f"log2 {x}"
            enrich_stats[new_x] = np.log2(enrich_stats[x])
            v = 0
        else:
            new_x = x

        if hue is None:
            n_clr = 1
        else:
            n_clr = len(set(enrich_stats[hue]))
        sns.scatterplot(data=enrich_stats, x=new_x, y=y, s=100, sizes=(20, 200),
                        palette=sns.color_palette("Set1", n_clr),
                        ax=curr_ax, hue=hue, size=size)

        curr_ax.axvline(x=v)
        # ax.plot([0.5], [0.5], transform=ax.transAxes, color='black')
        # plt.axis('square')
        curr_ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
                  borderaxespad=0)
        curr_ax.set_title(k)
        count+=1
    if f_save is not None:
        plt.tight_layout()
        plt.savefig(f_save, bbox_inches='tight', dpi=300)
    if to_close:
        plt.close()
    return

