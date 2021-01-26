import glob
from os.path import join
import pandas as pd
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import os
import pickle
from pandarallel import pandarallel
import numpy as np
from numpanpar import parallel_df as pardf
from Bio import SeqIO
import click


def create_af_tensors(concat_f, save_f, maxBP=16571):
    # Create a numpy tensor, of |cell|-|pos|-|nucs| shape, for both coverage along with average base quality
    # Cells are found in the concat file
    nucs = ["A", "C", "G", "T"]
    concat_dir = os.path.dirname(concat_f)
    # coverage_df = pd.read_csv(
    #     glob.glob(os.path.join(concat_dir, f"*.coverage.txt.gz"))[0])
    coverage_df = pd.read_csv(concat_f, header=None)
    cells = list(set(coverage_df[1].apply(lambda x: x.replace(".bam","")).values))

    coverage_tensor = np.zeros([len(cells), maxBP, len(nucs)])
    bq_tensor = np.zeros([len(cells), maxBP, len(nucs)])

    print(cells[:5])
    # Create key dictionary for cells, nucs,
    cells_d = {val:i for i, val in enumerate(cells)}
    nucs_d = {val:i for i, val in enumerate(nucs)}
    position_d = {val:i for i, val in enumerate(range(1,maxBP+1))}

    # Make a nucleotide-by-MT dataframe to get  overall allele frequencies for each position
    # nt_df.loc[:,:] = 0 # Initialize as 0s

    # # Dictionary of coverage per cell
    # nt_cov_dict = dict()
    for n in ["A", "C", "G", "T"]:
        print("Nucleotide: ", n)
        curr_f = glob.glob(os.path.join(concat_dir, f"*.{n}.txt.gz"))[0]
        df = pd.read_csv(curr_f, header=None)
        df.columns = ["Position", "Cell", "Coverage", "BQ"]
        df["Cell"] = df["Cell"].apply(lambda x: x.replace(".bam", ""))
        for ind, val in tqdm(df.iterrows()):
            coverage_tensor[cells_d[val["Cell"]], position_d[val["Position"]],nucs_d[n]] = val['Coverage']
            bq_tensor[cells_d[val["Cell"]], position_d[val["Position"]],nucs_d[n]] = val['BQ']

    pickle.dump([coverage_tensor, bq_tensor, cells_d, nucs_d, position_d], open(save_f, "wb"))
    return


def create_cell_distance():
    return


def fill_af_by_cell_loop(cell_df, coverage_dir, type="coverage"):
    """Cell series where indices are the positions
    Loads in the A/C/G/T/coverage files, then loops through the positions and fills in the alternative allele AF"""

    for cell, cell_series in tqdm(cell_df.iterrows()):
        cell_nucs = dict()
        for n in ["A", "C", "G", "T", "coverage"]:
            cell_nucs[n] = pd.read_csv(
                join(coverage_dir, "CB_" + cell + "." + n + ".txt"),
                header=None)
            if n == "coverage":
                cell_nucs[n].columns = ["Position", "Cell", "Coverage"]
            else:
                cell_nucs[n].columns = ["Position", "Cell", "Coverage",
                                        "BQ"]
            cell_nucs[n] = cell_nucs[n].set_index("Position",
                                                  drop=False)
            cell_nucs[n]["Cell"] = cell_nucs[n]["Cell"].apply(
                lambda x: x.replace(".bam", ""))
            cell_nucs[n]["Coverage"] = cell_nucs[n]["Coverage"].fillna(
                0)
        # Loop through each position and calculate AF frequency,
        # where the position is the last character of the index name
        for ind, val in tqdm(cell_series.iteritems()):
            pos = int(ind[:-1])
            n = ind[-1]
            # curr_nuc_df = nt_cov_dict[n]
            if pos in cell_nucs[n].index.values:
                if type == 'coverage':
                    cell_df.loc[cell, ind] = cell_nucs[n].loc[pos][
                                                 "Coverage"] / \
                                             cell_nucs["coverage"].loc[pos][
                                                 "Coverage"]
                elif type=="BQ":
                    cell_df.loc[cell, ind] = cell_nucs[n].loc[pos,"BQ"]
                else:
                    print("Type variable not recognized")

                if cell_df.loc[cell, ind] > 1:
                    print("Above 1", ind, pos)
            else:
                cell_df.loc[cell, ind] = 0

    return cell_df


@click.command(help="Create AF tensor (nCells-by-nPositions-by-nNucs) for coverage along with BQ tensor for average quality")
@click.argument('concat_f',  type=click.Path(exists=True))
@click.argument('af_f', type=click.Path(exists=False))
@click.option('-maxBP', type=click.INT, default=16571)
def main(concat_f, af_f, maxbp):
    "/data2/mito_lineage/data/processed/mttrace/CD34_mt_PolydT/mapq_30/scPileup_concat_200"
    create_af_tensors(concat_f, af_f, maxBP=maxbp)

if __name__ == '__main__':
    main()

