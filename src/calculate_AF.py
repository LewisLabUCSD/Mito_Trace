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


def calculate_af(coverage_dir, concat_dir, AF_F, ref_fasta, maxBP, topN=500, min_cells=100, min_reads=100, coverage_thresh=-1, het_thresh=0.05, min_het_cells=10):
    """

    :param coverage_dir:
    :param sc_coverage_f:
    :param AF_F:
    :param ref_fasta:
    :param maxBP:
    :param topN:
    :return:
           af_by_cell: DataFrame of allele-frequency of the top variant for each cell.
           af:  DataFrame with columns "Nucleotide", "AF", "Reference" and the rows being the number of MT positions.

    """
    nucs = ["A", "C", "G", "T"]
    mt_fasta_sequence = list(SeqIO.parse(ref_fasta, "fasta"))

    # pandarallel.initialize(nb_workers=32)

    # Create a numpy tensor, of |cell|-|pos|-|nucs| shape, for both coverage along with average base quality
    # Cells are found in the concat file
    coverage_df = pd.read_csv(glob.glob(os.path.join(concat_dir, f"*.coverage.txt.gz"))[0])
    #cells = coverage_df[1].values


    # Make a nucleotide-by-MT dataframe to get  overall allele frequencies for each position
    nt_df = pd.DataFrame(index=range(1, maxBP + 1), columns=nucs,
                         dtype=int)
    bq_df = pd.DataFrame(index=range(1, maxBP + 1), columns=nucs)

    # nt_df.loc[:,:] = 0 # Initialize as 0s

    # # Dictionary of coverage per cell
    # nt_cov_dict = dict()

    cells = dict()
    for n in ["A", "C", "G", "T"]:
        curr_f = glob.glob(os.path.join(concat_dir, f"*.{n}.txt.gz"))[0]
        df = pd.read_csv(curr_f, header=None)
        df.columns = ["Position", "Cell", "Coverage", "BQ"]
        df["Cell"] = df["Cell"].apply(lambda x: x.replace(".bam", ""))
        cells[n] = set(df["Cell"].values)
        pos = df.groupby("Position")
        nt_df.loc[:, n] = pos["Coverage"].agg("sum")
        bq_df.loc[:, n] = pos["BQ"].agg("mean")
        # nt_cov_dict[n] = df

    nt_df = nt_df.fillna(0)
    bq_df = bq_df.fillna(0)

    cell_names = set()
    for n in cells:
        cell_names=cell_names.union(cells[n])

    # af is Allele frequency dataframe to give the MT,
    # along with the variant allele frequency of the top variant and with reference allele
    af = pd.DataFrame(index=range(1, maxBP + 1),
                      columns=["Nucleotide", "AF", "Reference", "Alternative BQ"])
    no_alt_count = 0
    for ind, val in nt_df.iterrows():
        if ind > len(str(mt_fasta_sequence[0].seq)) - 1:
            break
        # Get the reference nucleotide from the fasta file
        curr_ref = str(mt_fasta_sequence[0].seq).upper()[ind - 1]

        if curr_ref not in nt_df.columns.values:
            continue
        alt_nuc = val.drop(curr_ref).idxmax() # Get the alternative allele
        if val[alt_nuc] == 0:
            no_alt_count += 1
            continue
        else:
            af.loc[ind, "Nucleotide"] = alt_nuc
            af.loc[ind, "AF"] = val[alt_nuc] / val.sum()
            af.loc[ind, "Reference"] = curr_ref
            # Store the BQ average value of the alternative allele
            af.loc[ind, "Alternative BQ"] = bq_df.loc[ind, alt_nuc]
    print('no alt count', no_alt_count)
    af = af.dropna()

    sc_coverage = pd.read_csv(glob.glob(join(concat_dir, "*.coverage.txt.gz"))[0],
                              header=None)
    sc_coverage.columns = ["Position", "Cell", "Coverage"]
    sc_coverage["Cell"] = sc_coverage["Cell"].apply(lambda x: x.replace(".bam",""))


    if topN == -1: ## Keep them all
        topCells = sc_coverage.groupby("Cell").sum()["Coverage"].sort_values(ascending=False)
    else:
        topCells = sc_coverage.groupby("Cell").sum()["Coverage"].sort_values(ascending=False)[:topN]
    cell_filter = topCells.index.values


    above_thresh = sc_coverage.groupby("Cell").sum()["Coverage"]>coverage_thresh
    above_thresh = above_thresh[above_thresh==True].index


    sc_coverage = sc_coverage.loc[sc_coverage["Cell"].isin(cell_filter)]

    # Only positions with minimum number of reads overall in the cell
    sc_coverage = sc_coverage.loc[sc_coverage["Cell"].isin(above_thresh.values)]

    # Only positions with minimum number of reads at the position
    pos_counts_filter = (sc_coverage.groupby("Position").apply(
        lambda x: (x['Coverage'] >= min_reads).sum())) >= min_cells
    pos_counts_filter = pos_counts_filter[pos_counts_filter==True].index.values

    print(f"Number of positions to keep : {len(pos_counts_filter)}")
    print(f"Number of cells to keep : {len(cell_filter)}")

    # Create the af_by_cell df, where index is the position concatenated with alternative allele
    af_by_cell = pd.DataFrame(
        index=af.apply(lambda x: (str(x.name) + x["Nucleotide"]),
                       axis=1).values, columns=cell_names)
    af_by_cell.columns = map(lambda x: x.replace("CB_", "").replace(".bam",""),
                             af_by_cell.columns)
    # Only keep filtered positions and cells
    af_by_cell = af_by_cell.loc[:, af_by_cell.columns.isin(cell_filter)]

    # Get the positions of variants that pass the filter
    af_by_cell_mask = af_by_cell.apply(lambda x: int(x.name[:-1]), axis=1)
    af_by_cell = af_by_cell.loc[
        af_by_cell_mask[af_by_cell_mask.isin(pos_counts_filter)].index]
    af_by_cell = af_by_cell.fillna(0)
    bq_by_cell = af_by_cell.copy()

    af_by_cell = pardf(af_by_cell.transpose(),
                            fill_af_by_cell_loop,
                            func_args=(coverage_dir,),
                            num_processes=32)
    af_by_cell = af_by_cell.loc[:,(af_by_cell>0).any(axis=0)]

    bq_by_cell = pardf(bq_by_cell.transpose(),
                            fill_af_by_cell_loop,
                            func_args=(coverage_dir,),
                            num_processes=32)
    bq_by_cell = bq_by_cell.loc[:,(bq_by_cell>0).any(axis=0)]

    het_filt = (af_by_cell > het_thresh).sum(axis=0) > min_het_cells
    het_filt = het_filt[het_filt].index
    af_by_cell = af_by_cell[het_filt]
    bq_by_cell = bq_by_cell[het_filt]

    if AF_F is not None:
        af_by_cell.to_csv(AF_F)
        bq_by_cell.to_csv(AF_F.replace('csv','')) + ".bq.csv"
    return af_by_cell, af, bq_df, nt_df


# def plot_af(af_by_cell_f, savedir):
#     af_by_cell = pd.read_csv(af_by_cell_f, index_col=0)
#     sns.clustermap(af_by_cell.fillna(0), row_cluster=False, col_cluster=False)
#     plt.savefig(os.path.join(savedir, "RAW_af_by_cell.png"))
#
#     top100_vars = af_by_cell.var().sort_values(ascending=False)[
#                   :100].index
#     sns.clustermap(np.sqrt(af_by_cell.loc[:, top100_vars]))
#     plt.savefig(os.path.join(savedir, "sqrt_top100Pos_af_by_cell.png"))
#     top100_vars = af_by_cell.var().sort_values(ascending=False)[:100].index
#     sns.clustermap(np.clip(np.sqrt(af_by_cell.loc[:,top100_vars]),a_min=0.05, a_max=0.4))
#     plt.savefig(os.path.join(savedir, "sqrt_clip_top100Pos_af_by_cell.png"))
#     return



@click.command(help="Create AF tensor (nCells-by-nPositions-by-nNucs) for coverage along with BQ tensor for average quality")
@click.argument('concat_f',  type=click.Path(exists=True))
@click.argument('af_f', type=click.Path(exists=False))
@click.option('-maxBP', type=click.INT, default=16571)
def main(concat_f, af_f, maxbp):
    "/data2/mito_lineage/data/processed/mttrace/CD34_mt_PolydT/mapq_30/scPileup_concat_200"
    create_af_tensors(concat_f, af_f, maxBP=maxbp)


def test():
    ref_fa = "/data2/mito_lineage/BWA-Primers-MT/MT_genome/MT.fasta"
    maxBP=16571
    concat_wt_dir="/data2/mito_lineage/data/processed/mttrace/CD34_mt_PolydT/mapq_30/scPileup_concat_200"
    coverage_wt_dir="/data2/mito_lineage/data/processed/mttrace/CD34_mt_PolydT/mapq_30/CD34_mt_PolydT_scPileup_200"
    calculate_af(coverage_wt_dir, concat_wt_dir, ref_fasta=ref_fa,
                 AF_F=None, maxBP=maxBP, topN=500, min_cells=100,
                 min_reads=10)

    return


if __name__ == '__main__':
    main()

