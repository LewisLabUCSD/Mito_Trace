import pandas as pd
import numpy as np
import pickle
import click
from sklearn.mixture import GaussianMixture

##### Filters
# def create_cell_filter(sc_coverage_f, filter_f, topN, coverage_count=0):
#     sc_coverage = pd.read_csv(sc_coverage_f,
#                               index_col=0)
#     top500 = sc_coverage.loc[
#         sc_coverage.sum(axis=1).sort_values(ascending=False)[
#         :topN].index]
#     cell_filter = top500.index.values
#     np.savetxt(fname=filter_f, X=cell_filter, fmt='s')
#     return

def create_cell_coverage():
    return


def create_barcode_filter(coverage_tensor, cells_d,bq_tensor,barcode_f):
    """
    Use the cells dictionary and the barcode names in barcodes_f to filter for the indicies of the tensors
    :param coverage_tensor:
    :param cells_d: {cell:index} dictionary for the tensors
    :param bq_tensor:
    :param barcode_f:
    :return:
    """
    CB_read_number = pickle.load(open(barcode_f, "rb"))
    cb_keys = list(CB_read_number.keys())

    cell_inds_to_keep = []
    cells_new_d = {}

    for k in cells_d:
        if k in cb_keys:
            #This will be the new index
            cells_new_d[k] = len(cell_inds_to_keep)

            cell_inds_to_keep.append(cells_d[k])

    coverage_tensor = coverage_tensor[cell_inds_to_keep]
    bq_tensor= bq_tensor[cell_inds_to_keep]


    return coverage_tensor, bq_tensor, cells_new_d


def create_position_counts_filter(coverage_tensor, bq_tensor, pos_d, min_cells=100, min_reads=100):
    #pos_filter = (coverage_tensor.sum(axis=2) >= min_reads).sum(axis=0) >= min_cells
    pos_filter = np.where((coverage_tensor.sum(axis=2) >= min_reads).sum(
        axis=0) >= min_cells)[0]
    coverage_tensor = coverage_tensor[:, pos_filter,:]
    bq_tensor = bq_tensor[:, pos_filter, :]

    new_pos_d = {}
    for p in pos_d:
        if pos_d[p] in pos_filter: #if in filters inds, add as the new index
            new_pos_d[p] = len(new_pos_d)

    return coverage_tensor, bq_tensor, new_pos_d


def create_top_positions_filter(coverage_tensor,bq_tensor, pos_dict, topPositions=-1):
    coverage_pos = coverage_tensor.sum(axis=2).sum(axis=0)
    pos_inds = coverage_pos.argsort()[::-1]
    if topPositions  > 0:
        pos_inds = pos_inds[:topPositions]

    #pos = {c:pos_dict[c] for c in pos_inds}
    new_pos_d = {}
    for c in pos_dict:
        if pos_dict[c] in pos_inds: #if in filters inds, add as the new index
            new_pos_d[c] = len(new_pos_d)


    coverage_tensor = coverage_tensor[:,pos_inds, :]
    bq_tensor = bq_tensor[:, pos_inds, :]
    return coverage_tensor, bq_tensor, new_pos_d


def create_top_cells_filter(coverage_tensor, bq_tensor, cell_dict, topCells=-1, by='coverage'):
    #| cell | - | pos | - | nucs |
    coverage_cells = coverage_tensor.sum(axis=1).sum(axis=1)
    cell_inds = coverage_cells.argsort()[::-1]
    if topCells  > 0:
        cell_inds = cell_inds[:topCells]

    #cells_d = {c:cell_dict[c] for c in cell_inds}
    new_cell_d = {}
    for c in cell_dict:
        if cell_dict[c] in cell_inds: #if in filters inds, add as the new index
            new_cell_d[c] = len(new_cell_d)

    coverage_tensor = coverage_tensor[cell_inds, :, :]
    bq_tensor = bq_tensor[cell_inds, :, :]

    return coverage_tensor, bq_tensor, new_cell_d


def bq_filter(coverage_tensor, bq_tensor, pos_d, type='Gaussian'):
    bq_vals = bq_tensor.flatten()
    gm = GaussianMixture(n_components=3)
    gm.fit(bq_vals)
    vals=gm.get_params()

def af_distance():
    return


##### Till here
@click.command(help="Filter AF and BQ tensor (nCells-by-nPositions-by-nNucs)")
@click.argument('af_tensor_f',  type=click.Path(exists=True))
@click.argument('barcode_f', type=click.Path(exists=False))
@click.argument('out_af_filter_f', type=click.Path(exists=False))
@click.argument('minCells', type=click.INT, default=-1)
@click.argument('minPos', type=click.INT, default=-1)
@click.argument('topCells', type=click.INT, default=-1)
@click.argument('topPos', type=click.INT, default=-1)
@click.argument('cell_mt_coverage', type=click.INT, default=-1)
def create_af_cell(af_tensor_f, barcode_f, out_af_filter_f, mincells, minpos, topcells, toppos, cell_mt_coverage):
    # Load the AF data
    [coverage_tensor, bq_tensor, cells_d, nucs_d, position_d] = pickle.load(open(af_tensor_f,"rb"))

    # Create the filters step by step

    # Barcode filter
    coverage_tensor, bq_tensor, cells_d = create_barcode_filter(coverage_tensor, cells_d, bq_tensor, barcode_f)

    # Top cells coverage
    coverage_tensor, bq_tensor, cells_d = create_top_cells_filter(coverage_tensor, bq_tensor, cells_d, topCells=topcells)

    # Top positions
    coverage_tensor, bq_tensor, position_d = create_top_positions_filter(coverage_tensor, bq_tensor, position_d, topPositions=toppos)

    # Positions with minimum number of reads by number of cells
    coverage_tensor, bq_tensor, position_d = create_position_counts_filter(coverage_tensor, bq_tensor, position_d, min_cells=mincells, min_reads=minpos)

    pickle.dump([coverage_tensor, bq_tensor, cells_d, position_d],open(out_af_filter_f, 'wb'))
    return


if __name__ == '__main__':
    create_af_cell()