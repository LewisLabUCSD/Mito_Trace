import pandas as pd
import glob
import os
from tqdm import tqdm
import pickle
from numpanpar import parallel_df as pardf
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys
import click
import matplotlib as mpl
mpl.use('agg')


def fill_df_coverage(df, pileup_dir, is_par=False):
    """ Reads each cell file and stores the coverage for each cell"""
    not_done = []
    df = df.astype(int)
    for ind in (df.index):
        f = glob.glob(os.path.join(pileup_dir,"CB_" + ind + ".coverage.txt"))
        f_rev = glob.glob(os.path.join(pileup_dir,"CB_" + ind + ".coverage.minus.txt"))

        if len(f) == 0 and len(f_rev) == 0:
            not_done.append(ind)
        else:
            curr = pd.read_csv(f[0], header=None)
            for _, val in curr.iterrows():
                df.loc[ind, val[0]] = val[2]
            if len(f_rev) > 0:  # Reverse strand
                curr_rev = pd.read_csv(f_rev[0], header=None)
                for _, val in curr_rev.iterrows():
                    df.loc[ind, val[0]] += val[2]

    print(f"Total {df.shape[0] - len(not_done)}")
    print(f"Number of missing files: {len(not_done)}")
    return df



# @click.command(help="Create MT by position")
# @click.option('--maxbp', default=16571, help='The length of the MT genome')
# @click.argument('barcode_p',  type=click.Path(exists=True))
# @click.argument('pileup_dir',  type=click.Path(exists=True))
# @click.argument('save_f', type=click.STRING)
def sc_mt_coverage(barcode_p, concat_coverage_f, save_f, maxBP):
    """ Converts from sparse to wide, and uses the Cell Barcodes

    :param barcode_p: The cellbarcode information file. Counts how many reads per barcode
    :param concat_coverage_f: The concatenated, long, format of all cells
    :param save_f:
    :param maxBP: This is the cutoff number of basepairs needed to be covered to be added to the filtered matrix
    :return:
    """
    print(barcode_p, concat_coverage_f, save_f, maxBP)
    CB_read_number = pickle.load(open(barcode_p, "rb"))

    df = pd.read_csv(concat_coverage_f, header=None)
    # Only use the + strand here
    sc_coverage = (df.pivot(index=1, columns=0, values=2)).fillna(0)
    sc_coverage.index = sc_coverage.index.str.replace(".bam", "")
    sc_coverage = sc_coverage.loc[sc_coverage.index.isin(list(CB_read_number.keys()))]
    sc_coverage.to_csv(save_f)
    return


def plot_sc_mt(sc_coverage_f, savefig_f, top_n=0):
    sc_coverage = pd.read_csv(sc_coverage_f, index_col=0)
    # top500 = sc_coverage.loc[
    #     sc_coverage.sum(axis=1).sort_values(ascending=False)[
    #     :top_n].index]
    # # sns.clustermap(top500, col_cluster=False)
    # sns.clustermap(top500, col_cluster=False, vmax=500)

    if not ((top_n == 0) or (top_n==-1)):
        sc_coverage = sc_coverage.loc[
            sc_coverage.sum(axis=1).sort_values(ascending=False)[
            :min(top_n, sc_coverage.shape[0])].index]
    else:
        top_n = sc_coverage.shape[0]
    # else:
    #     sc_coverage = sc_coverage.loc[
    #         sc_coverage.sum(axis=1).sort_values(ascending=False)]

    sns.heatmap(sc_coverage) #, row_cluster=False, col_cluster=False)
    plt.title(f"Number of reads at each MT position in {top_n} cells")
    #g.ax_heatmap.set_title(f"Number of reads at each MT position by the top {top_n} covered cells")
    #g.ax_heatmap.set_yticks([])
    #g.ax_heatmap.set_xlabel("MT position")
    plt.xlabel("MT position")
    #plt.legend(title="Coverage (Number of reads)")
    plt.ylabel("Cell")
    #g.ax_heatmap.set_ylabel("Cell")
    plt.savefig(savefig_f)
    plt.savefig(savefig_f.replace(".png", ".svg"))
    log2_sc_coverage = np.log2(sc_coverage + 1)
    if not top_n == 0:
        log2_sc_coverage = log2_sc_coverage.loc[
            log2_sc_coverage.sum(axis=1).sort_values(ascending=False)[
            :top_n].index]
    g = sns.clustermap(log2_sc_coverage, col_cluster=False)
    g.ax_heatmap.set_title(f"Log2 number of reads at each MT position in {top_n} cells")
    g.ax_heatmap.set_yticks([])
    g.ax_heatmap.set_xlabel("MT position")
    g.ax_heatmap.set_ylabel("Cell")
    plt.savefig(savefig_f.replace(".png","")+".log2.png")
    plt.savefig(savefig_f.replace(".png", ".svg"))

    return


def plot_percent_coverage_cells(sc_coverage_f, savefig_f,
                                x_cov=(1, 2, 5, 10, 30, 50, 100, 200, 300),
                                cells_cov = (1, 10, 100, 500, 1000,5000,10000)):
    sc_coverage = pd.read_csv(sc_coverage_f, index_col=0)
    # x_cov = [1, 2, 5, 10, 30, 50, 100, 200, 1000, 2000, 3000]
    # cells_cov = [1, 10, 100, 500]

    pct_cov_df = pd.DataFrame(
        columns=["Number of MT positions", "Minimum number of cells",
                 "Minimum number of reads"])

    for x in x_cov:
        curr_count = (sc_coverage >= x).sum(axis=0)
        for c in cells_cov:
            pct_cov_df = pd.concat((pct_cov_df, pd.DataFrame(
                {"Number of MT positions": (curr_count >= c).sum(),
                 "Minimum number of cells": c,
                 "Minimum number of reads": x}, index=[len(
                    pct_cov_df)])))  # pct_cov.append((sc_coverage > x).sum(axis=0).sum())
    f = plt.figure()
    sns.barplot(data=pct_cov_df, x="Minimum number of reads",
                y="Number of MT positions",
                hue="Minimum number of cells")
    plt.savefig(savefig_f)
    plt.savefig(savefig_f.replace(".png",".svg"))
    return



# @click.command(help="Create MT by position")
# @click.option('--maxbp', default=16571, help='The length of the MT genome')
# @click.argument('barcode_p',  type=click.Path(exists=True))
# @click.argument('pileup_dir',  type=click.Path(exists=True))
# @click.argument('save_f', type=click.STRING)
def main(func, *args, **kwargs):
    print('args', args)
    if func == "sc_mt":
        barcode_p, concat_coverage_f, save_f, maxBP = args
        maxBP = int(maxBP)
        #print(barcode_p, pileup_dir, save_f, maxBP)
        sc_mt_coverage(barcode_p, concat_coverage_f, save_f, maxBP)#, n_cpu=)

    elif func == "plot":
        sc_coverage_f, save_f_pos_cov = args
        #sc_coverage_f, save_f_pos_heat, save_f_pos_cov = args
        if "x_cov" in kwargs:
            x_cov = kwargs["x_cov"]
        else:
            x_cov = [1, 2, 5, 10, 30, 50, 100, 200, 1000, 2000, 3000]
        if "cells_cov" in kwargs:
            cells_cov = kwargs["cells_cov"]
        else:
            cells_cov = [1, 10, 100, 500, 1000,5000,10000]
        plot_percent_coverage_cells(sc_coverage_f, save_f_pos_cov , x_cov, cells_cov)
        #plot_sc_mt(sc_coverage_f, save_f_pos_heat, top_n=1000)
    return


if __name__ == "__main__":
    #sc_mt_coverage()
    main(sys.argv[1], *sys.argv[2:])

