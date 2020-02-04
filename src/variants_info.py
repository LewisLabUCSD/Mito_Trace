import glob
from os.path import join
import pandas as pd
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import click
import os
import pickle
from pandarallel import pandarallel
import numpy as np

# To remove variants whose inferred heteroplasmy may reflect sequencing errors, we examined the distribution of per-base, per-allele base-quality scores, noting a clear pattern of high quality and low-quality variants (Figure S1C). To determine high quality variants, we fit a mixture of three Gaussian distributions (Figure S1C, labeled by different colors), and filtered such that only alleles that had > 99% probability of belonging to the blue (largest mean BQ) Gaussian were retained. This conservative empirical threshold for a BQ cutoff was determined to be 23.8 based on this mixture model approach (Figure S1C, vertical dotted line). As one poorly quantified position allele would affect the estimates for all other alleles at the specific position, we filtered positions
# that contained one or more alleles with a BQ less than the empirical threshold unless the allele had a non-significant (i.e., less than 1 in 500) effect on heteroplasmy. In total, we called 44 high-quality variants across our TF1 (sub-)clones (Figure S1D) that were present at a minimum of 2.5% heteroplasmy in at least one sample.

def create_filtered_variants(coverage_folder, ref_fasta, BQ_file, base_quality_threshold=0):
    """
    Function that gets variants to be used in constructing distance matrix

    For each variant at each position, see if
    1) Number of reads is cumulatively higher than number of reads
    2) %TO DO% BQ threshold

    :return:

    """
    BQ = pd.read_csv(coverage_folder, sep="\t")

    return


def par_variant_fill(bq_init_val):
    #for f in BQ_init.index.values:
    f = bq_init_val.name
    curr = pd.read_csv(f, header=None)
    curr.columns = ["Position", "Cell", "Coverage", "BQ"]
    for _, val in curr.iterrows():
        bq_init_val[val["Position"]] = val["BQ"]
        #bq_init_val.at[f, val["Position"]] = val["BQ"]
    return bq_init_val


def par_plot_variant_quality(coverage_folder, f_save, maxBP = 16571, f_save_fig=None):
    pandarallel.initialize(nb_workers=32)
    A_files = glob.glob(join(coverage_folder, "*.A.txt"))
    C_files = glob.glob(join(coverage_folder, "*.C.txt"))
    G_files = glob.glob(join(coverage_folder, "*.G.txt"))
    T_files = glob.glob(join(coverage_folder, "*.T.txt"))
    nuc_files = {"A": A_files,"C":C_files,"G":G_files,"T":T_files}

    #cells = list(map(lambda x: os.path.basename(x.replace(".T.txt","")),T_files))
    BQ = pd.DataFrame(index=range(1,maxBP+1), columns=["A","C","G","T"],dtype=int)

    nuc_df_dict = {}
    for n in tqdm(nuc_files):
        BQ_init = pd.DataFrame(index=nuc_files[n],
                               columns=list(range(1, maxBP + 1)), dtype=int)

        BQ_init = BQ_init.iloc[:32,:]
        nuc_df_dict[n] = BQ_init.parallel_apply(par_variant_fill,axis=1)

        BQ[n] = BQ_init.transpose().median(axis=1, skipna=True)

    # Save files
    pickle.dump(obj=nuc_df_dict, file=open(f_save + ".p","wb"))
    BQ.to_csv(f_save + ".median.tsv", sep="\t")
    #BQ_cells.to_csv(f_save.replace(".tsv","")+".cells.tsv", sep="\t")

    f, ax = plt.subplots()
    bq_values = BQ.values.flatten()
    bq_values = bq_values[~np.isnan(bq_values)]
    ax.hist(bq_values)
    if f_save_fig is None:
        plt.savefig(f_save+".png")
    else:
        plt.savefig(f_save_fig)
    return




def create_fullsum_matrix(allele_folder, f_save,groups=None):
    """
    Function to prepare data for downstream processing from mito-genotyping github page

    Output:
        n-cells-by-4 tsv file with header as "Cell", "group", "barcodeSequences", "mtCoverage
    :return :
    """

    return


def create_distance(variants):
    return


# def cluster_cells(cell_position_nuc_f,variants_list):
#     cell_position_nuc= pickle.load((open(cell_position_nuc_f,"rb")))
#     for n in nuc:
#         nuc[n]
#     pd.read_csv()
#     return


# @click.command(help="Will plot BQ distribution")
# @click.option('--maxBP', default=16571, help='The length of the MT genome')
# @click.option('--f_save_fig', default=None, help='The folder we see')
# @click.argument('coverage_folder',  type=click.Path(exists=True)) #help='file to append the command',
# @click.argument('f_save', type=click.STRING)#coverage_folder, f_save, maxBP = 16571, f_save_fig=None
def main(coverage_folder, f_save, maxBP, f_save_fig):
    # plot_variant_quality(coverage_folder, f_save, maxBP=16571,
    #                      f_save_fig=None)
    par_plot_variant_quality(coverage_folder, f_save, maxBP=16571,
                         f_save_fig=None)
    return


if __name__ == "__main__":
    coverage_folder = "/data2/mito_lineage/data/processed/A/scPileup_concat_200"
    f_save = "/data2/mito_lineage/data/processed/A/DownstreamProcessing/BQ"
    maxBP = 16571
    main(coverage_folder=coverage_folder, f_save=f_save, maxBP=maxBP, f_save_fig=None)

