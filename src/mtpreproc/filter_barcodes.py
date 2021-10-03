import click
import pickle
import pandas as pd
from os.path import join, basename


def filter_barcodes_from_CB(CB_read_number, cellr_bc_f):
    """
    Function that filters a dictionary where the keys are cells from a text file of cell barcodes, outputted from cellranger
    :param CB_read_number: dictionary where the keys are the cell IDs
    :param cellr_bc_f: Each line is a cell ID
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


@click.command(help="Filter cell barcodes based on cellranger list")
@click.argument('barcode_f',  type=click.Path(exists=True))
@click.argument('cellr_f', type=click.STRING)
@click.argument('use_cellr', type=click.BOOL)
@click.argument('out_f', type=click.Path())
def main(barcode_f, cellr_f, use_cellr, out_f):
    """ Filters barcodes. If a list is given, use that, otherwise load the whole set from the MT bam file
    Args:
        barcode_f:
        cellr_f:
        use_cellr:
        out_f:
    """
    if cellr_f == "" or use_cellr is False:
        cellr_f = None

    [_, CB_read_number, _, _,
     _, _] = pickle.load(open(barcode_f, "rb"))
    if cellr_f is not None:
        CB_read_number = filter_barcodes_from_CB(CB_read_number, cellr_f)
    pickle.dump(CB_read_number,file=open(out_f, "wb"))
    with open(out_f.replace(".p", ".txt"), "w") as f:
        f.write("\n".join([f"CB:Z:{x}" for x in list(CB_read_number.keys())]))

    return



if __name__== '__main__':
    main()
