import pickle
import numpy as np

def test_create_barcode_filter(coverage_tensor, cells_d,bq_tensor,barcode_f):
    """Use the cells dictionary and the barcode names in barcodes_f to filter
    for the indicies of the tensors :param coverage_tensor: :param cells_d:
    {cell:index} dictionary for the tensors :param bq_tensor: :param barcode_f:
    :return:

    Args:
        coverage_tensor:
        cells_d:
        bq_tensor:
        barcode_f:
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


def test_create_position_counts_filter(coverage_tensor, bq_tensor, pos_d, min_cells=100, min_reads=100):
    #pos_filter = (coverage_tensor.sum(axis=2) >= min_reads).sum(axis=0) >= min_cells
    """
    Args:
        coverage_tensor:
        bq_tensor:
        pos_d:
        min_cells:
        min_reads:
    """
    pos_filter = np.where((coverage_tensor.sum(axis=2) >= min_reads).sum(
        axis=0) >= min_cells)[0]
    coverage_tensor = coverage_tensor[:, pos_filter,:]
    bq_tensor = bq_tensor[:, pos_filter, :]

    new_pos_d = {}
    for p in pos_d:
        if pos_d[p] in pos_filter: #if in filters inds, add as the new index
            new_pos_d[p] = len(new_pos_d)

    return coverage_tensor, bq_tensor, new_pos_d


def test_create_top_positions_filter(coverage_tensor,bq_tensor, pos_dict, topPositions=-1):
    """
    Args:
        coverage_tensor:
        bq_tensor:
        pos_dict:
        topPositions:
    """
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


def test_create_top_cells_filter(coverage_tensor, bq_tensor, cell_dict, topCells=-1, by='coverage'):
    #| cell | - | pos | - | nucs |
    """
    Args:
        coverage_tensor:
        bq_tensor:
        cell_dict:
        topCells:
        by:
    """
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
    """
    Args:
        coverage_tensor:
        bq_tensor:
        pos_d:
        type:
    """
    bq_vals = bq_tensor.flatten()
    gm = GaussianMixture(n_components=3)
    gm.fit(bq_vals)
    vals=gm.get_params()

def af_distance():
    return