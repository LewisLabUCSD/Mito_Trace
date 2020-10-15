import pandas as pd
from os.path import basename
import numpy as np
import seaborn as sns
from mplh.fig_utils import helper_save
import matplotlib.pyplot as plt


def scatter_vaf(sample_df, save_f=None, hue=None):
    """
    Paired scatterplot across all the samples for heteroplasmy
    :param sample_df:
    :param save_f:
    :return:
    """
    sns.pairplot(sample_df, hue=hue, vars=list(filter(lambda x: x!= hue, sample_df.columns.values)))
    plt.title("Variant Allele Frequency")
    helper_save(save_f)
    #plt.savefig("Variant Allele Frequency")
    return


def run_compare_vaf(sample_list, inds_f, save_f):
    #sample_dict = dict()
    overlap_variants = pd.read_csv(inds_f, header=None)
    samples = list(map(lambda x: basename(x), sample_list))

    sample_df = pd.DataFrame(index=overlap_variants, columns=samples, dtype=np.float)
    for ind, f in enumerate(sample_list):

        sample_df.loc[overlap_variants , samples[ind]] = pd.read_csv(f)["AF"]

    scatter_vaf(sample_df, save_f)
    return


def plot_top_hets(af_by_cell, hets, sample_col=None, f_save=None):
    f, ax = plt.subplots(nrows=len(hets))
    for ind, h in enumerate(hets):
        sns.violinplot(af_by_cell[h], hue=sample_col, inner = "stick", ax=ax[ind])
    helper_save(f_save)
    return
