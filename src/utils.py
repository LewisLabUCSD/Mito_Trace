import pandas as pd


def filter_barcodes_from_CB(CB_read_number, cellr_bc_f):
    """
    Function that filters a dictionary where the keys are cells from a text file of cell barcodes, outputted from cellranger
    :param CB_read_number:
    :param cellr_bc_f:
    :return:
    """
    if cellr_bc_f is None:
        return CB_read_number
    filt_bc = pd.read_csv(cellr_bc_f, sep='\t', header=None)[0].values
    count = 0
    CB_read_number_filt = {}
    for x in filt_bc:
        if x not in set(CB_read_number.keys()):
            count += 1
            print(x)
        else:
            CB_read_number_filt[x] = CB_read_number[x]

    print(f"Number of missing barcodes in CB but present in filtered: {count}")
    print(f"Number of cells after using cellranger cell filter: {len(CB_read_number_filt)}")
    return CB_read_number_filt
