import glob
import os
import pandas as pd
from tqdm import tqdm
import numpy as np
import pysam
import time
from collections import defaultdict
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import subprocess
from os.path import join

# FUNC_PATH = "/home/isshamie/software/mito-genotyping/exampleProcessing/"
# BAM = "../data/raw/14515X1.MT.bam"
# BARCODES = "../data/raw/Human2_002_Barcodes.txt"
# NUM_CORES = 32
# IN = "../data/raw/scBam"
# OUT_DIR = "../data/raw/scPileup/"
# BARCODE_INFO = "../data/scBam/barcode_data.p"
# maxBP = 16571
#
# if not os.path.exists(OUT_DIR):
#     os.mkdir(OUT_DIR)
#
# files = glob.glob(IN + "/*")
# print(len(files))

## Load CB's, get ones that have over average read coverage of 1x
def plot_coverage(barcode_data, f_save, maxBP=16571):
    CB_read_number = pickle.load(open(barcode_data, "rb"))

    CB_read_MT = dict()
    for i in CB_read_number:
        if CB_read_number[i] * 100 >= maxBP:
            CB_read_MT[i] = CB_read_number[i] * 100

    plt.figure()
    plt.hist(np.array(list(CB_read_MT.values())) / maxBP, log=True)
    plt.xlabel("X Coverage of Mitochondrion")
    plt.ylabel("Number of cells")
    plt.title(f"LOG Number of nts per MT length per cell\nN cells with > 1x read coverage ={len(CB_read_MT)}")

    plt.savefig(f_save) #"../figures/CB_coverage_hist.png")
    return


def plot_coverage_barcode_file(barcode_data, barcode_fs, f_save, maxBP=16571):
    """
    Plots the coverage of MT by cells that are in a barcode file(s)
    :param barcode_data: pickled barcode info, outputted from split_by_CB
    :param barcode_fs: list of barcode files to be read in
    :param f_save: Figure name to save
    :param maxBP: Length of the MT
    :return:
    """
    if type(barcode_fs) == list:
        barcode_text = pd.read_csv(barcode_fs[0], sep="\t", index_col=0)
        for i in barcode_fs[1:]:
            barcode_text = pd.concat((barcode_text, pd.read_csv(i, sep="\t", index_col=0)))
    else:
        barcode_text = pd.read_csv(barcode_fs, sep="\t", index_col=0)

    CB_read_number = pickle.load(open(barcode_data, "rb"))

    CB_read_MT = dict()
    cb = 0
    for i in barcode_text["Barcode"].values:
        i = i.replace(".", "-")
        # print(i)
        if i in CB_read_number:
            curr = barcode_text[
                barcode_text["Barcode"] == i.replace("-", ".")]
            # print(barcode_text[barcode_text["Barcode"] == i.replace("-",".")])
            # print(i)

            print(curr["Sample.Name"].values[0])
            if curr["Sample.Name"].values[0] == "Human2_002":
                print("Here")
            cb += 1
            if CB_read_number[i]*100 >= maxBP:
                CB_read_MT[i] = CB_read_number[i]*100

    plt.figure()
    plt.hist(np.array(list(CB_read_MT.values())) / maxBP, log=True)
    plt.xlabel("X Coverage of Mitochondrion")
    plt.ylabel("Number of cells")
    plt.title(
        f"LOG Number of nts per MT length per cell\nN cells with > 1x read coverage ={len(CB_read_MT)}")

    plt.savefig(f_save)
    return


def plot_coverage_from_scPileup(scPileup_f, barcodes_f, f_save, cov_thresh=1, maxBP=16571):
    scPileup  = pd.read_csv(scPileup_f, header=None)
    scPileup.columns = ["Position", "Cell", "Coverage"]
    cell_coverage = scPileup.groupby("Cell").sum()["Coverage"]
    cell_coverage_mt = cell_coverage/maxBP
    cell_coverage_mt.index = list(map(lambda x: x.replace('.bam',''),cell_coverage_mt.index.values))


    CB = list(pickle.load(open(barcodes_f, "rb")).keys())
    print(cell_coverage_mt.index)
    print('before cb filter', len(cell_coverage_mt.index))
    cell_coverage_mt = cell_coverage_mt.loc[cell_coverage_mt.index.isin(CB)]
    print('after cb filter', len(cell_coverage_mt.index))

    print('cov thresh', cov_thresh)
    cell_coverage_mt = cell_coverage_mt[cell_coverage_mt>cov_thresh]
    print('after cov filter', len(cell_coverage_mt.index))




    plt.figure()
    plt.hist(cell_coverage_mt, log=True)
    plt.xlabel("Coverage of Mitochondrion Genome")
    plt.ylabel("Number of cells")
    #plt.title(f"MT coverage\nN cells with > {cov_thresh}x read coverage ={len(cell_coverage_mt)}")
    plt.title(f"MT coverage per cell")
    plt.savefig(f_save)

    f = plt.figure()
    sns.distplot(cell_coverage_mt)
    plt.savefig(f_save.replace('.png','.dist.png'))

    return


if __name__ == "__main__":
    plot_coverage_from_scPileup(sys.argv[1], sys.argv[2], sys.argv[3], maxBP=16571)
    #plot_coverage(sys.argv[1], sys.argv[2], maxBP=16571)
