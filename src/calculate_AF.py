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
from numpanpar import parallel_df as pardf
from Bio import SeqIO


def fill_af_by_cell(cell_series, coverage_folder):
    """Cell series where indices are the positions
    Loads in the A/C/G/T/coverage files, then loops through the positions and fills them in"""
    cell = cell_series.name
    print('cell',cell)
    cell_nucs = dict()
    for n in ["A","C","G","T","coverage"]:
        cell_nucs[n] = pd.read_csv(join(coverage_folder,"CB_" + cell + "." + n + ".txt"), header=None)

        if n == "coverage":
                cell_nucs[n].columns = ["Position", "Cell", "Coverage"]
        else:
            cell_nucs[n].columns = ["Position", "Cell", "Coverage", "BQ"]
        cell_nucs[n] = cell_nucs[n].set_index("Position", drop=False)
        cell_nucs[n]["Cell"] = cell_nucs[n]["Cell"].apply(lambda x: x.replace(".bam",""))
        cell_nucs[n]["Coverage"] =  cell_nucs[n]["Coverage"].fillna(0)
    # Loop through each position and calculate AF frequency
    for ind, val in tqdm(cell_series.iteritems()):
        pos = int(ind[:-1])
        n = ind[-1]
        #curr_nuc_df = nt_cov_dict[n]
        if pos in cell_nucs[n].index.values:
            cell_series.loc[ind] = cell_nucs[n].loc[pos]["Coverage"] / cell_nucs["coverage"].loc[pos]["Coverage"]
            if cell_series.loc[ind] > 1:
                print("Above 1", ind, pos)
        else:
            cell_series.loc[ind] = 0

    return cell_series


def fill_af_by_cell_loop(cell_df, coverage_folder):
    """Cell series where indices are the positions
    Loads in the A/C/G/T/coverage files, then loops through the positions and fills them in"""

    for cell, cell_series in tqdm(cell_df.iterrows()):
        cell_nucs = dict()
        for n in ["A", "C", "G", "T", "coverage"]:
            cell_nucs[n] = pd.read_csv(
                join(coverage_folder, "CB_" + cell + "." + n + ".txt"),
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
        # Loop through each position and calculate AF frequency
        for ind, val in tqdm(cell_series.iteritems()):
            pos = int(ind[:-1])
            n = ind[-1]
            # curr_nuc_df = nt_cov_dict[n]
            if pos in cell_nucs[n].index.values:
                cell_df.loc[cell, ind] = cell_nucs[n].loc[pos][
                                             "Coverage"] / \
                                         cell_nucs["coverage"].loc[pos][
                                             "Coverage"]
                if cell_df.loc[cell, ind] > 1:
                    print("Above 1", ind, pos)
            else:
                cell_df.loc[cell, ind] = 0

    return cell_df


def create_cell_filter(sc_coverage_f, filter_f, topN, coverage_count=0):
    sc_coverage = pd.read_csv(sc_coverage_f,
                              index_col=0)
    top500 = sc_coverage.loc[
        sc_coverage.sum(axis=1).sort_values(ascending=False)[
        :topN].index]
    cell_filter = top500.index.values
    np.savetxt(fname=filter_f, X=cell_filter, fmt='s')
    return


def create_position_counts_filter(sc_coverage_f, filter_f, min_cells=100, min_reads=100):
    sc_coverage = pd.read_csv(sc_coverage_f,
                              index_col=0)
    #min_cells, min_reads = 100, 100
    counts_filter = (sc_coverage >= min_reads).sum(axis=0) >= min_cells
    counts_filter = counts_filter[counts_filter == True].index.values
    np.savetxt(fname=filter_f, X=counts_filter, fmt='s')
    return counts_filter


def calculate_af(coverage_folder, sc_coverage_f, AF_F,  ref_fasta, maxBP, topN=500):
    nucs = ["A", "C", "G", "T"]
    records = list(SeqIO.parse(ref_fasta, "fasta"))
    # pandarallel.initialize(nb_workers=32)
    A_files = glob.glob(join(coverage_folder, "*.A.txt"))
    C_files = glob.glob(join(coverage_folder, "*.C.txt"))
    G_files = glob.glob(join(coverage_folder, "*.G.txt"))
    T_files = glob.glob(join(coverage_folder, "*.T.txt"))
    nuc_files = {"A": A_files, "C": C_files, "G": G_files, "T": T_files}

    cells = list(
        map(lambda x: os.path.basename(x.replace(".T.txt", "")),
            T_files))

    nt_df = pd.DataFrame(index=range(1, maxBP + 1), columns=nucs,
                         dtype=int)
    # nt_df.loc[:,:] = 0 # Initialize as 0s
    af = pd.DataFrame(index=range(1, maxBP + 1),
                      columns=["Nucleotide", "AF", "Reference"])
    nt_cov_dict = dict()

    for n in ["A", "C", "G", "T"]:
        curr_f = glob.glob(os.path.join(coverage_folder, f"*.{n}.txt.gz"))[0]
        df = pd.read_csv(curr_f, header=None)
        df.columns = ["Position", "Cell", "Coverage", "BQ"]
        df["Cell"] = df["Cell"].apply(lambda x: x.replace(".bam", ""))
        pos = df.groupby("Position")
        nt_df.loc[:, n] = pos["Coverage"].agg("sum")
        nt_cov_dict[n] = df

    nt_cov_dict["coverage"] = pd.read_csv(
        glob.glob(os.path.join(coverage_folder, "*.coverage.txt.gz"))[
            0], header=None)
    nt_cov_dict["coverage"].columns = ["Position", "Cell", "Coverage"]
    nt_cov_dict["coverage"]["Cell"] = nt_cov_dict["coverage"][
        "Cell"].apply(lambda x: x.replace(".bam", ""))

    nt_df = nt_df.fillna(0)

    for ind, val in nt_df.iterrows():
        if ind > len(str(records[0].seq)) - 1:
            break
        curr_ref = str(records[0].seq).upper()[ind - 1]

        if curr_ref not in nt_df.columns.values:
            continue
        alt_nuc = val.drop(curr_ref).idxmax()
        af.loc[ind, "Nucleotide"] = alt_nuc
        af.loc[ind, "AF"] = val[alt_nuc] / val.sum()
        af.loc[ind, "Reference"] = curr_ref

    af = af.dropna()
    sc_coverage = pd.read_csv(sc_coverage_f,
                              index_col=0)
    top500 = sc_coverage.loc[
        sc_coverage.sum(axis=1).sort_values(ascending=False)[
        :topN].index]
    cell_filter = top500.index.values

    min_cells, min_reads = 100, 100
    counts_filter = (top500 >= min_reads).sum(axis=0) >= min_cells
    counts_filter = counts_filter[counts_filter == True].index.values


    print(f"Number of positions to keep : {len(counts_filter)}")
    print(f"Number of cells to keep : {len(cell_filter)}")

    af_by_cell = pd.DataFrame(
        index=af.apply(lambda x: (str(x.name) + x["Nucleotide"]),
                       axis=1).values, columns=cells)
    af_by_cell.columns = map(lambda x: x.replace("CB_", ""),
                             af_by_cell.columns)
    # Only keep filtered positions and cells
    af_by_cell = af_by_cell.loc[:, af_by_cell.columns.isin(cell_filter)]

    af_by_cell_mask = af_by_cell.apply(lambda x: x.name[:-1], axis=1)
    af_by_cell = af_by_cell.loc[
        af_by_cell_mask[af_by_cell_mask.isin(counts_filter)].index]
    af_by_cell = af_by_cell.fillna(0)

    af_by_cell = pardf(af_by_cell.transpose(),
                            fill_af_by_cell_loop,
                            func_args=(coverage_folder,),
                            num_processes=32)

    af_by_cell.to_csv(AF_F)
    return


def plot_af(af_by_cell_f, savedir):
    af_by_cell = pd.read_csv(af_by_cell_f, index_col=0)
    sns.clustermap(af_by_cell.fillna(0), row_cluster=False, col_cluster=False)
    plt.savefig(os.path.join(savedir, "RAW_af_by_cell.png"))

    top100_vars = af_by_cell.var().sort_values(ascending=False)[
                  :100].index
    sns.clustermap(np.sqrt(af_by_cell.loc[:, top100_vars]))
    plt.savefig(os.path.join(savedir, "sqrt_top100Pos_af_by_cell.png"))
    top100_vars = af_by_cell.var().sort_values(ascending=False)[:100].index
    sns.clustermap(np.clip(np.sqrt(af_by_cell.loc[:,top100_vars]),a_min=0.05, a_max=0.4))
    plt.savefig(os.path.join(savedir, "sqrt_clip_top100Pos_af_by_cell.png"))
    return
