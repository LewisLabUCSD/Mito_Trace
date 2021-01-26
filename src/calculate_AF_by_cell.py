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


def fill_af_by_cell_loop(cell_df, coverage_dir, type="coverage"):
    """ Cell series where indices are the positions and columns are variant info.

    Loads in the A/C/G/T/coverage files, then loops through the
    positions and fills in the alternative allele AF.
    """

    for cell, cell_series in cell_df.iterrows():
        cell_nucs = dict()
        for n in ["A", "C", "G", "T", "coverage"]:
            if not os.path.exists(join(coverage_dir, "CB_" + cell + "." + n + ".txt")): #it only has the rev strand
                print(f"{n} not found for {cell}")
                continue
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
        for ind, val in (cell_series.iteritems()):
            pos = int(ind[:-1])
            n = ind[-1]
            # curr_nuc_df = nt_cov_dict[n]
            if n not in cell_nucs or "coverage" not in cell_nucs:
                continue
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


def filter_positions_cells(sc_coverage, topN=500, min_cells=100, min_reads=100, coverage_thresh=-1):
    """ Filters cells and counts based on depth

    :param sc_coverage:
    :param topN:
    :param min_cells:
    :param min_reads:
    :param coverage_thresh:
    :return:
    """
    if topN == -1 or topN == 0: ## Keep them all
        print("topN", topN)
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
    return pos_counts_filter, cell_filter


def filter_allele_cells(af_by_cell, het_thresh=0, min_het_cells=0):
    """ Filter the variants based on counts or AF thresholds.

    :param af_by_cell:
    :type af_by_cell: pd.DataFrame or np.array of floats or ints
    :param het_thresh:
    :param min_het_cells:
    :return:
    """
    het_filt = (af_by_cell > het_thresh).sum(axis=0) > min_het_cells
    het_filt = het_filt[het_filt].index
    print(f'Positions that pass het filter: {len(het_filt)}')
    return het_filt


def get_nt_bq(concat_dir, maxBP):
    """ For each position, get the counts and average quality for each nt.

    :param concat_dir:
    :param maxBP:
    :return:
    """
    nucs = ["A", "C", "G", "T"]
    # Make a |nucleotide|-by-MT dataframe to get  overall allele frequencies for each position
    nt_df = pd.DataFrame(index=range(1, maxBP + 1), columns=nucs,
                         dtype=int)
    bq_df = pd.DataFrame(index=range(1, maxBP + 1), columns=nucs)

    # nt_df.loc[:,:] = 0 # Initialize as 0s

    # # Dictionary of coverage per cell
    # nt_cov_dict = dict()

    cells = dict()
    for n in ["A", "C", "G", "T"]:
        curr_f = glob.glob(os.path.join(concat_dir, f"*.{n}.strands.txt.gz"))[0]
        df = pd.read_csv(curr_f)

        # Only take the Forward reads since it is the one to compare against the reference
        if len(df.columns)>4:
            df = df.iloc[:, :4]
        df.columns = ["Position", "Cell", "Coverage", "BQ"]
        df["Cell"] = df["Cell"].apply(lambda x: x.replace(".bam", ""))
        cells[n] = set(df["Cell"].values)
        pos = df.groupby("Position")
        nt_df.loc[:, n] = pos["Coverage"].agg("sum")
        bq_df.loc[:, n] = pos["BQ"].agg("mean")
        # nt_cov_dict[n] = df

    nt_df = nt_df.fillna(0)
    bq_df = bq_df.fillna(0)

    return nt_df, bq_df, cells


def filter_and_save_nt_pileups(concat_dir, cell_inds, nt_pos, out_d, name):
    """ Filters the NT pileup and bq matrices and saves as pileups.

    :param concat_dir:
    :param cell_inds:
    :param nt_pos:
    :param out_d:
    :param name:
    :return:
    """
    for n in ["A", "C", "G", "T"]:
        curr_f = glob.glob(os.path.join(concat_dir, f"*.{n}.strands.txt.gz"))[0]
        df = pd.read_csv(curr_f)
        # Only take the Forward reads since it is the one to compare against the reference
        if len(df.columns)>4:
            df.columns = ["Position", "Cell", "+ Coverage", "+ BQ", "- Coverage", "- BQ"]
            #df = df.iloc[:, :4]
        else:
            df.columns = ["Position", "Cell", "Coverage", "BQ"]
        nt_df = nt_df[nt_df["Cell"].isin(cell_inds)]
        nt_df = nt_df[nt_df["Position"].isin(nt_pos)]
        df["Cell"] = df["Cell"].apply(lambda x: x.replace(".bam", ""))

        curr_out = join(out_d, f"{name}.{n}.txt")
        df.to_csv(curr_out, header=None, index=False)
    return


def load_sc_pileup(sc_f):
    """ Loads and sets the columns of the coverage pileup.
    :param sc_f:
    :return:
    """
    sc_coverage = pd.read_csv(glob.glob(sc_f)[0])
    if len(sc_coverage)>3:
        sc_coverage = sc_coverage.iloc[:, :3]
    sc_coverage.columns = ["Position", "Cell", "Coverage"]
    sc_coverage["Cell"] = sc_coverage["Cell"].apply(lambda x: x.replace(".bam",""))
    return sc_coverage


def extract_af(nt_df, bq_df, mt_fasta_sequence, maxBP):
    """ Creates the position alleles and their quality

    :param nt_df:
    :param bq_df:
    :param mt_fasta_sequence:
    :param maxBP:
    :return:
    """
    af = pd.DataFrame(index=range(1, maxBP + 1),
                      columns=["Nucleotide", "AF", "Reference",
                               "Alternative BQ"])
    no_alt_count = 0
    for ind, val in nt_df.iterrows():
        if ind > len(str(mt_fasta_sequence[0].seq)) - 1:
            break
        # Get the reference nucleotide from the fasta file
        curr_ref = str(mt_fasta_sequence[0].seq).upper()[ind - 1]

        if curr_ref not in nt_df.columns.values:
            continue
        alt_nuc = val.drop(
            curr_ref).idxmax()  # Get the alternative allele
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

    return af


def calculate_af(coverage_dir, concat_dir, AF_F, ref_fasta, maxBP,
                 topN=500, min_cells=100, min_reads=100,
                 coverage_thresh=-1, het_thresh=0, min_het_cells=0,
                 num_proc=8, bq_thresh=0, het_count_thresh=0, name=""):
    """ Creates a filtered allele frequency matrix and saves in different ways.

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

    nt_df, bq_df, cells = get_nt_bq(concat_dir, maxBP)
    # pandarallel.initialize(nb_workers=32)

    cell_names = set()
    for n in cells:
        cell_names=cell_names.union(cells[n])

    # Load the MTPos-by-cell matrix
    sc_coverage = load_sc_pileup(join(concat_dir, "*.coverage.strands.txt.gz"))


    # FILTER: top cells by coverage, minimum coverage, minimum number of cells and certain number of reads
    pos_counts_filter, cell_filter = filter_positions_cells(sc_coverage, topN=topN, min_cells=min_cells,
                           min_reads=min_reads, coverage_thresh=coverage_thresh)

    if len(cell_filter) == 0 or len(pos_counts_filter) == 0:
        print("Nothing passes the filter.")
        if AF_F is not None:
            os.system(f'touch {AF_F}')
            os.system(f'touch {AF_F.replace("csv", "") + ".bq.csv"}')
            #bq_by_cell.to_csv(AF_F.replace('csv', '')) + ".bq.csv"
            return None

    # af is the AF Pos-by-[AF, BQ] df
    mt_fasta_sequence = list(SeqIO.parse(ref_fasta, "fasta"))
    af = extract_af(nt_df, bq_df, mt_fasta_sequence, maxBP)

    ## Filter here
    ###############
    # Create the af_by_cell df, where index is the position concatenated with alternative allele (e.g. 13C) and column is cell barcode
    af_by_cell = pd.DataFrame(
        index=af.apply(lambda x: (str(x.name) + x["Nucleotide"]),
                       axis=1).values, columns=cell_names)
    af_by_cell.columns = list(map(lambda x: x.replace("CB_", "").replace(".bam",""),
                             af_by_cell.columns))
    # Only keep filtered positions and cells
    af_by_cell = af_by_cell.loc[:, af_by_cell.columns.isin(cell_filter)]

    # Get the positions of variants that pass the filter
    af_by_cell_mask = af_by_cell.apply(lambda x: int(x.name[:-1]), axis=1)
    af_by_cell = af_by_cell.loc[
        af_by_cell_mask[af_by_cell_mask.isin(pos_counts_filter)].index]
    af_by_cell = af_by_cell.fillna(0)
    #af_by_cell_mask = af_by_cell_mask.loc[af_by_cell.index]
    #Make the bq_by_cell, the same as af_by_cell except average BQ instead of summed AF coverage
    ###############

    af_by_cell = pardf(af_by_cell.transpose(),
                            fill_af_by_cell_loop,
                            func_args=(coverage_dir,),
                            num_processes=num_proc)
    af_by_cell = af_by_cell.loc[:,(af_by_cell>0).any(axis=0)]

    bq_by_cell = af_by_cell.copy() #Placeholder for shape
    bq_by_cell = pardf(bq_by_cell.transpose(),
                            fill_af_by_cell_loop,
                            func_args=(coverage_dir,),
                            num_processes=num_proc)

    bq_by_cell = bq_by_cell.loc[:,(bq_by_cell>0).any(axis=0)]

    # FILTER by minimum cells w heterozygous AF minimum
    het_filt = filter_allele_cells(af_by_cell, het_thresh, min_het_cells=min_het_cells)
    af_by_cell = af_by_cell[het_filt]
    bq_by_cell = bq_by_cell[het_filt]


    # FILTER by minimum cells w heterozygous count minimum
    het_filt = filter_allele_cells(af_by_cell, het_count_thresh, min_het_cells=min_het_cells)
    af_by_cell = af_by_cell[het_filt]
    bq_by_cell = bq_by_cell[het_filt]


    # FILTER by average BQ after filtering.
    bq_average = bq_by_cell.mean(axis=0)
    bq_avg_filt = bq_by_cell[bq_average >= bq_thresh].columns
    af_by_cell = af_by_cell[bq_avg_filt]


    # Save AF and BQ matrix, but also
    if AF_F is not None:
        af_by_cell.to_csv(AF_F)
        bq_by_cell.to_csv(AF_F.replace('.csv','') + ".bq.csv")

        # Save filtered nucleotide pileups for mgatk usage
        filter_and_save_nt_pileups(concat_dir, af_by_cell.columns, pos_counts_filter, os.path.dirname(AF_F), name)

        # Save depth and allele frequency as sparse matrices
        #af_by_cell.melt()

    return af_by_cell, bq_by_cell, af, bq_df, nt_df



@click.command(help="Create AF by cell (nCells-by-nPositions-by-nNucs) for coverage along with BQ tensor for average quality")
@click.argument('coverage_dir',  type=click.STRING)
@click.argument('concat_dir', type=click.STRING)
@click.argument('af_f', type=click.STRING)
@click.argument('min_cells', type=click.INT)
@click.argument('min_reads', type=click.INT)
@click.argument('topn', type=click.INT)
@click.argument('coverage_thresh',type= click.INT)
@click.argument('het_thresh', type=click.FLOAT)
@click.argument('min_het_cells', type=click.INT)
@click.argument('maxbp', type=click.INT)
@click.argument('ref_fa', type=click.Path(exists=True))
def main(coverage_dir, concat_dir, af_f, min_cells, min_reads,topn, coverage_thresh, het_thresh,min_het_cells,maxbp,ref_fa):
    #ref_fa = "/data2/mito_lineage/BWA-Primers-MT/MT_genome/MT.fasta"
    #maxbp = 16571
    calculate_af(coverage_dir, concat_dir, af_f, ref_fa, maxBP=maxbp,
                 topN=topn, min_cells=min_cells,
                 min_reads=min_reads, coverage_thresh=coverage_thresh,
                 het_thresh=het_thresh, min_het_cells=min_het_cells)


# @click.command(help="Create AF by cell (nCells-by-nPositions-by-nNucs) for coverage along with BQ tensor for average quality")
# @click.argument('coveragedir',  type=click.STRING)
# @click.argument('concat_dir', type=click.STRING)
# @click.argument('af_f', type=click.STRING)
# @click.argument('min_cells', type=click.INT)
# @click.argument('min_reads', type=click.INT)
# @click.argument('topn', type=click.INT)
# @click.argument('coveragethresh',type= click.INT)
# def main(coveragedir, concat_dir, af_f, min_cells, min_reads,topn, coveragethresh ):
#     print(coveragedir)
#     print(af_f)
#     print(topn)
#     print(min_cells)
#     print(coveragethresh)
#     return


if __name__ == '__main__':
    main()

