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
    [CR_read_number, CB_read_number, BC_read_number, barcodes,
     corrected_barcodes, barcode_pairs] = pickle.load(
        open(barcode_data, "rb"))


    CB_read_MT = dict()
    for i in CB_read_number:
        if CB_read_number[i] * 100 >= maxBP:
            CB_read_MT[i] = CB_read_number[i] * 100

    plt.figure()
    plt.hist(np.array(list(CB_read_MT.values())) / maxBP, log=True)
    plt.xlabel("X Coverage of Mitochondrion")
    plt.ylabel("Number of cells")
    plt.title(
        f"LOG Number of nts per MT length per cell\nN cells with > 1x read coverage ={len(CB_read_MT)}")

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

    [CR_read_number, CB_read_number, BC_read_number, barcodes,
     corrected_barcodes, barcode_pairs] = pickle.load(
        open(barcode_data, "rb"))

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


if __name__ == "__main__":
    plot_coverage(sys.argv[1], sys.argv[2], maxBP=16571)
