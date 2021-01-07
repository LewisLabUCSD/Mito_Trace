from os.path import join
import os
from numpy import random
import numpy as np
import pandas as pd


def run_batch(indirs, outdir, num_reads_total):
    num_reads_each = len(indirs)/num_reads_total

    return


def subsample_sparse_matrices(indirs, outdir, cell_subsample=0.1,
                              equalCellNumbers=False, num_cells_total=1000,
                              variant_subsample=0):
    """ Subsamples the cellSNP output matrices and combines into a new folder, with metadata saying which cell comes from which sample.

    :param indirs:
    :param outdir:
    :param cell_subsample:
    :param equalCellNumbers:
    :param num_cells_total:
    :param variant_subsample:
    :return:
    """
    for i in indirs:
        ad_f = join(i, "cellSNP.tag.AD.mtx")
        dp_f = join(i, "cellSNP.tag.DP.mtx")
        oth_f = join(i, "cellSNP.tag.OTH.mtx")

        ad = pd.read_csv(ad_f, comment="%", header=None, sep="\t")
        ad.columns = ["Variant", "Cell", "integer"]
        ad = ad.iloc[1:]
        total_ad = ad["integer"].sum()

        #
        num_cells = cell_subsample*total_ad
    return


def merge_vcf(vcf_files, out_vcf):
    cmd = f"""bcftools merge {vcf_files} - Oz - o {out_vcf}"""
    print(cmd)
    os.system(cmd)
    return

# Subsampling on bam files:
# Pros: More raw
# Cons: Will subsampling on bam miss the paired reads? Can that be connected

#Subsampling on mtx files:
#Pros: Already takes care of the bam processing, and takes care of batch
# for each separately.
