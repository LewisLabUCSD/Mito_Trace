# try:
#     from .utils import utils
# except ImportError:
#     from src.utils import utils

from src.utils import utils
import click
import pickle


@click.command(help="Filter cell barcodes based on cellranger list")
@click.argument('barcode_f',  type=click.Path(exists=True))
@click.argument('cellr_f', type=click.STRING)
@click.argument('use_cellr', type=click.BOOL)
@click.argument('out_f', type=click.Path())
def main(barcode_f, cellr_f, use_cellr, out_f):
    # bam_f = sys.argv[1]
    # out_f = sys.argv[2]
    """
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
        CB_read_number = utils.filter_barcodes_from_CB(CB_read_number, cellr_f)

    pickle.dump(CB_read_number,file=open(out_f, "wb"))

    return



if __name__== '__main__':
    main()
